"""Driving the virtualPCR JAR.

One JAR invocation per (primer pair, target file). That is deliberate: given
several primers at once the JAR pairs *any* forward with *any* reverse and
reports chimeric amplicons across unrelated pairs. Its FRpairs mode avoids
that but stops reporting amplicon sequences and Ta. Running one pair at a time
keeps attribution unambiguous and the report rich; a full mitogenome set costs
a few seconds per pair, and pairs run in parallel.
"""
from __future__ import annotations

import os
import shutil
import struct
import subprocess
import tempfile
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

from .model import PairResult, PrimerPair, RunConfig
from .parser import fasta_headers, parse_report, repair_sequences
from .primers import write_pair_fasta

# Java class file major version -> feature release (45 == Java 1.1)
_CLASS_MAJOR_BASE = 44


class JarError(RuntimeError):
    pass


@dataclass(frozen=True)
class JavaInfo:
    executable: str
    feature: int          # e.g. 24
    raw: str


def _probe_java(exe: str) -> JavaInfo | None:
    try:
        p = subprocess.run([exe, "-version"], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return None
    raw = (p.stderr or p.stdout).strip()
    feature = 0
    for tok in raw.replace('"', " ").split():
        head = tok.split(".")[0]
        if head.isdigit():
            feature = int(head)
            break
    return JavaInfo(exe, feature, raw.splitlines()[0] if raw else "")


def find_java(explicit: str | None = None) -> JavaInfo:
    exe = explicit or os.environ.get("JAVA_HOME_BIN") or shutil.which("java")
    if not exe:
        raise JarError("java not found on PATH; install a JDK or set the path explicitly")
    info = _probe_java(exe)
    if info is None:
        raise JarError(f"could not run {exe}")
    return info


def discover_javas() -> list[JavaInfo]:
    """Every runnable JRE we can find, newest first."""
    seen: dict[str, JavaInfo] = {}
    cands: list[str] = []
    if (w := shutil.which("java")):
        cands.append(w)
    for root in (r"C:\Program Files\Java", r"C:\Program Files\Eclipse Adoptium",
                 "/usr/lib/jvm", "/Library/Java/JavaVirtualMachines"):
        p = Path(root)
        if not p.is_dir():
            continue
        for jdk in p.iterdir():
            for rel in ("bin/java.exe", "bin/java", "Contents/Home/bin/java"):
                exe = jdk / rel
                if exe.exists():
                    cands.append(str(exe))
    for exe in cands:
        if exe in seen:
            continue
        if (info := _probe_java(exe)):
            seen[exe] = info
    return sorted(seen.values(), key=lambda j: j.feature, reverse=True)


def select_java(jar: str | Path, explicit: str | None = None) -> JavaInfo:
    """A JRE new enough to run this JAR, preferring an explicit choice."""
    if explicit:
        java = find_java(explicit)
        check_compatibility(jar, java)
        return java
    need = jar_class_version(jar)
    java = find_java()
    if java.feature >= need:
        return java
    for cand in discover_javas():
        if cand.feature >= need:
            return cand
    raise JarError(
        f"{Path(jar).name} needs Java {need}+ but only Java {java.feature} was found. "
        f"Install a newer JDK, or use a JAR built for Java {java.feature}."
    )


def jar_class_version(jar: str | Path) -> int:
    """Feature release the JAR was compiled for, from its main class."""
    with zipfile.ZipFile(jar) as z:
        name = next((n for n in z.namelist() if n.endswith("virtualpcr.class")), None)
        if name is None:
            raise JarError(f"{jar} does not contain virtualpcr.class")
        head = z.read(name)[:8]
    if head[:4] != b"\xca\xfe\xba\xbe":
        raise JarError("not a Java class file")
    (major,) = struct.unpack(">H", head[6:8])
    return major - _CLASS_MAJOR_BASE


def check_compatibility(jar: str | Path, java: JavaInfo) -> None:
    need = jar_class_version(jar)
    if java.feature and java.feature < need:
        raise JarError(
            f"{Path(jar).name} needs Java {need}+ but '{java.executable}' is Java "
            f"{java.feature}. Install a newer JDK or use a JAR built for Java {java.feature}."
        )


def _link_or_copy(src: Path, dst: Path) -> None:
    """Give the JAR its own target path without duplicating a large FASTA."""
    try:
        os.link(src, dst)
        return
    except (OSError, NotImplementedError):
        pass
    try:
        os.symlink(src, dst)
        return
    except (OSError, NotImplementedError, AttributeError):
        pass
    shutil.copy2(src, dst)


@dataclass
class RunTask:
    pair: PrimerPair
    target: Path


def _safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in name) or "pair"


# The JAR announces failures two ways: an uncaught exception on stderr, or a
# plain 'Failed to ...' line on stdout while still exiting 0. Echoing the raw
# tail of stdout instead would bury the cause under the config it just parsed.
_ERROR_PREFIXES = ("Error", "Failed", "Incorrect", "Warning", "Skipped")


def _diagnose(stderr: str | None, stdout: str | None) -> str:
    for text in (stderr or "", stdout or ""):
        for line in text.splitlines():
            s = line.strip()
            if "Exception" in s or s.startswith(_ERROR_PREFIXES):
                return s
    tail = (stderr or stdout or "").strip()
    return tail[-200:] if tail else "no output"


def run_one(
    task: RunTask,
    jar: Path,
    java: JavaInfo,
    config: RunConfig,
    heap: str | None = None,
    keep_report: Path | None = None,
    headers: list[str] | None = None,
    repair: bool = True,
) -> PairResult:
    """Run the JAR for one pair against one FASTA and parse the report."""
    # The temp dir must sit on the target's volume so the hard link succeeds and
    # the JAR (older builds ignore output_path) writes its .out beside the link.
    parent = task.target.parent
    tmp = Path(tempfile.mkdtemp(prefix=".vpcr_", dir=parent))
    try:
        stem = _safe(task.pair.name)
        target_link = tmp / f"{stem}{task.target.suffix or '.fasta'}"
        _link_or_copy(task.target, target_link)

        primer_file = tmp / f"{stem}.primers.txt"
        write_pair_fasta(task.pair, primer_file)

        out_file = tmp / f"{stem}.out"
        cfg_file = tmp / f"{stem}.config"
        cfg_file.write_text(
            config.to_config_text(str(target_link), str(primer_file), str(out_file)),
            encoding="utf-8",
        )

        cmd = [java.executable, "-jar"]
        if heap:
            cmd[1:1] = [f"-Xmx{heap}"]
        cmd += [str(jar), str(cfg_file)]
        proc = subprocess.run(cmd, capture_output=True, text=True, cwd=tmp)

        # Older builds ignore output_path and write <target>.out next to it.
        report = out_file if out_file.exists() else target_link.with_suffix(
            target_link.suffix + ".out"
        )
        if not report.exists():
            raise JarError(
                f"no report produced for pair {task.pair.name!r} "
                f"(exit {proc.returncode}): {_diagnose(proc.stderr, proc.stdout)}"
            )

        result = parse_report(report, task.pair, target_file=str(task.target),
                              all_headers=headers)
        if repair:
            repair_sequences(result, task.target)
        if keep_report:
            keep_report.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(report, keep_report)
        return result
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def run_batch(
    pairs: Iterable[PrimerPair],
    targets: Iterable[str | Path],
    jar: str | Path,
    config: RunConfig | None = None,
    *,
    java_path: str | None = None,
    workers: int = 4,
    heap: str | None = None,
    reports_dir: str | Path | None = None,
    repair: bool = True,
    progress: Callable[[int, int, str], None] | None = None,
    should_cancel: Callable[[], bool] | None = None,
) -> list[PairResult]:
    """Run every pair against every target file, in parallel.

    `progress(done, total, label)` is called from worker threads.
    """
    # Resolved: workers run the JAR with cwd set to a temp dir, so any relative
    # path given by the caller would be interpreted against the wrong directory.
    jar = Path(jar).resolve()
    if not jar.exists():
        raise JarError(f"JAR not found: {jar}")
    java = select_java(jar, java_path)
    config = config or RunConfig()

    targets = [Path(t).resolve() for t in targets]
    for t in targets:
        if not t.exists():
            raise JarError(f"target file not found: {t}")
    pairs = list(pairs)
    if not pairs:
        raise JarError("no primer pairs given")

    header_cache = {t: fasta_headers(t) for t in targets}
    tasks = [RunTask(p, t) for t in targets for p in pairs]
    total = len(tasks)
    reports_dir = Path(reports_dir) if reports_dir else None

    # Keyed by task index so results can be restored to submission order:
    # completion order depends on thread scheduling, and callers derive folder
    # names and row order from this list.
    collected: dict[int, PairResult] = {}
    done = 0
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        futures = {}
        for i, t in enumerate(tasks):
            if should_cancel and should_cancel():
                break
            keep = (reports_dir / t.target.stem / f"{_safe(t.pair.name)}.out") if reports_dir else None
            futures[pool.submit(
                run_one, t, jar, java, config, heap, keep, header_cache[t.target], repair
            )] = (i, t)

        for fut in as_completed(futures):
            i, t = futures[fut]
            done += 1
            if should_cancel and should_cancel():
                fut.cancel()
                continue
            try:
                collected[i] = fut.result()
            except Exception as e:  # one bad pair must not kill the batch
                if progress:
                    progress(done, total, f"ERROR {t.pair.name}: {e}")
                continue
            if progress:
                progress(done, total, f"{t.pair.name} / {t.target.name}")
    return [collected[i] for i in sorted(collected)]
