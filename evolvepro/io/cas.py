from __future__ import annotations

import os
import shutil
import stat
import tempfile
from dataclasses import dataclass
from pathlib import Path

from evolvepro.io.hashing import sha256_file
from evolvepro.io.schemas import FileHash


@dataclass(frozen=True)
class Snapshot:
    path: Path
    file_hash: FileHash


def _snapshot_path(cas_root: Path, digest: str) -> Path:
    return cas_root / digest[:2] / digest[2:4] / digest


def _make_read_only(path: Path) -> None:
    path.chmod(stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)


def snapshot_file(source: Path, cas_root: Path) -> Snapshot:
    source = source.expanduser()
    digest = sha256_file(source)
    file_hash = FileHash(path=str(source), sha256=digest, size=source.stat().st_size)
    destination = _snapshot_path(cas_root, digest)

    if destination.exists():
        if sha256_file(destination) != digest:
            raise ValueError(f"CAS hash mismatch for existing snapshot: {destination}")
        _make_read_only(destination)
        return Snapshot(path=destination, file_hash=file_hash)

    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "wb", dir=destination.parent, delete=False
    ) as handle:
        tmp_path = Path(handle.name)
        with source.open("rb") as source_handle:
            shutil.copyfileobj(source_handle, handle)
        handle.flush()
        os.fsync(handle.fileno())

    _make_read_only(tmp_path)
    tmp_path.replace(destination)
    return Snapshot(path=destination, file_hash=file_hash)
