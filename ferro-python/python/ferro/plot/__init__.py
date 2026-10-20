"""Publication figures from ferro's csv products: ``gr``, ``sq``, ``angle``, ``msd``, ``net``.

One module per figure, each with a ``plot(inputs, outdir=None, suffix=None, ...)`` that writes
``<name>.pdf``, ``<name>.png`` and ``<name>_data.csv`` and returns their paths; ``ferro-py plot``
is the command line over them. ``common`` is the shared layer (style, reading the csv,
``labels.csv``, the data export).

Nothing is imported here: ``import ferro`` must not pull in matplotlib, which is an optional
dependency (``pip install 'ferro[plot]'``, plus LaTeX).
"""
