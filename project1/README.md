### Introduction

This is my project one slideshow setup. The presentation itself should contain enough information to recreate FEA in Ansys or another software.


### Installation

First, install LaTeX. You can find more information about how to install it [here](https://www.latex-project.org/get/).


This slideshow uses the Gotham theme. You can either install and use this theme, or simply comment out the line about the theme and pick another one in `presentation.tex`. You can find a guide to install it [here](https://gitlab.com/RomainNOEL/beamertheme-gotham#how-to-install-).


### Compiling

Compile using your preferred GUI, or use the command line. You might need to change your config settings to produce a PDF. You will need to be in this (`project1`) directory to do this.
```bash
latexmk presentation.tex
```
You should now see a file called `presentation.pdf`. If it has a different file extension (that isn't `.tex`), change your LaTeX engine settings. If you open it and see (??) instead of an equation number, simply compile the document again and it will update.
