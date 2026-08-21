# Presentations

Advisor-facing slide decks and their reproducible build scripts live here.

- `HCC_RNAseq_downstream_advisor_20260821.pptx` is the six-slide Chinese summary of the RNA-seq downstream analysis.
- `HCC_RNAseq_downstream_advisor_20260821.html` is the animated, offline-capable browser version of the same six-slide narrative.
- `HCC_RNAseq_downstream_advisor_20260821_standalone.html` embeds every image and can be downloaded and opened directly without a web server or SSH tunnel.
- `build_rnaseq_deck.py` regenerates the deck and presentation-specific chart assets from the tracked result tables.
- `make_standalone_html.py` regenerates the single-file HTML after the source HTML or images change.
- `assets/` contains compact figures used by the deck; they may be regenerated at any time.

Run from the project root:

```bash
.venv/bin/python presentations/build_rnaseq_deck.py
python3 presentations/make_standalone_html.py
```

For a remote workspace, the simplest option is to download
`HCC_RNAseq_downstream_advisor_20260821_standalone.html` and double-click it on
the local computer. It does not require `localhost`, an HTTP server, or an SSH
tunnel.

Open the animated version directly by double-clicking the HTML file. For the most reliable full-screen experience, serve the project locally:

```bash
cd /home/cx/DrugRep
python3 -m http.server 8000
```

If the browser is running on the same machine, open:

`http://127.0.0.1:8000/presentations/HCC_RNAseq_downstream_advisor_20260821.html`

If the HTTP server is running on a remote machine, create an SSH tunnel from a
terminal on the local computer. The SSH port must be supplied with `-p` (do not
append it to the hostname):

```bash
ssh -p 22222 -N -L 18000:127.0.0.1:8000 cx@120.76.203.127
```

Then open:

`http://127.0.0.1:18000/presentations/HCC_RNAseq_downstream_advisor_20260821.html`

If a local SSH proxy configuration incorrectly redirects the connection to a
closed proxy port, bypass the SSH configuration for this connection:

```bash
ssh -F /dev/null -o ProxyCommand=none -o ProxyJump=none \
  -p 22222 -N -L 18000:127.0.0.1:8000 cx@120.76.203.127
```

Keep both the remote HTTP-server terminal and the local SSH-tunnel terminal
open while presenting.
