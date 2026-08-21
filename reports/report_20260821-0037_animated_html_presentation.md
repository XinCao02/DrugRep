# Animated HTML advisor presentation

Created a local, offline-capable HTML version of the six-slide RNA-seq advisor deck:

[`HCC_RNAseq_downstream_advisor_20260821.html`](../presentations/HCC_RNAseq_downstream_advisor_20260821.html)

It preserves the PowerPoint narrative, restrained ivory/navy/teal/gold styling, Chinese-first language, and `ZZH` naming convention. Added interactions are limited to presentation needs: staged entrance animations, animated DEG counts, workflow progression, slide transitions, progress display, keyboard and mouse navigation, touch swiping, and browser full-screen mode. Reduced-motion system settings disable nonessential motion.

The page has no CDN or network dependency. It loads the five validated figures from `presentations/assets/` and can be opened directly or through a local HTTP server.

Final Chromium checks at 1600×900 confirmed one active slide at a time, zero horizontal or vertical overflow on all six slides, successful loading of every image, working arrow-key navigation, and no browser runtime errors.

Recommended launch:

```bash
cd /home/cx/DrugRep
python3 -m http.server 8000
```

Open `http://localhost:8000/presentations/HCC_RNAseq_downstream_advisor_20260821.html`. Use left/right arrows or Space to navigate and `F` for full screen.
