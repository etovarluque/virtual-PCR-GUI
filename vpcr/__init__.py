"""In silico PCR over FASTA targets, driving the virtualPCR JAR."""
from pathlib import Path

__version__ = "1.0.0"

APP_NAME = "virtualPCR Studio"

# Folder holding Run vPCR.bat, README.md and jar/. Used to locate the JAR.
APP_ROOT = Path(__file__).resolve().parent.parent
JAR_DIR = APP_ROOT / "jar"
