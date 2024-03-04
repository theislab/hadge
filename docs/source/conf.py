project = "hadge"
copyright = "2023, Fabiola Curion, Xichen Wu, Lukas Heumos"
author = "Fabiola Curion, Xichen Wu, Lukas Heumos"
release = "1.0.0"

extensions = ["myst_parser", "nbsphinx"]

source_suffix = {
    ".rst": "restructuredtext",
    ".txt": "markdown",
    ".md": "markdown",
}

templates_path = ["_templates"]

exclude_patterns = ["_build", "**.ipynb_checkpoints"]

html_theme = "furo"
html_static_path = ["_static"]

nbsphinx_execute = "never"
