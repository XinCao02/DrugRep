#!/usr/bin/env python3
"""Embed presentation images so the HTML opens without HTTP or SSH."""

from __future__ import annotations

import base64
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "HCC_RNAseq_downstream_advisor_20260821.html"
OUTPUT = ROOT / "HCC_RNAseq_downstream_advisor_20260821_standalone.html"
IMAGE_PATTERN = re.compile(r'src="(assets/[^"?]+\.png)"')


def embed_image(match: re.Match[str]) -> str:
    relative_path = match.group(1)
    image_path = ROOT / relative_path
    if not image_path.is_file():
        raise FileNotFoundError(f"Missing presentation image: {image_path}")
    payload = base64.b64encode(image_path.read_bytes()).decode("ascii")
    return f'src="data:image/png;base64,{payload}"'


def main() -> None:
    source_html = SOURCE.read_text(encoding="utf-8")
    embedded_html, replacement_count = IMAGE_PATTERN.subn(embed_image, source_html)
    if replacement_count != 5:
        raise RuntimeError(f"Expected 5 embedded images, found {replacement_count}")
    OUTPUT.write_text(embedded_html, encoding="utf-8")
    print(f"Created {OUTPUT} with {replacement_count} embedded images")


if __name__ == "__main__":
    main()
