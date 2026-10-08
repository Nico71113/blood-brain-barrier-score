#!/bin/zsh
set -e
cd "${0:A:h}"
if [[ ! -x .venv/bin/python ]]; then
  python3 -m venv .venv
  .venv/bin/python -m pip install -r requirements.txt
fi
print 'Open http://127.0.0.1:8000. Press Control-C to stop.'
exec .venv/bin/python app.py
