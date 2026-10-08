"""A local-only web app. Run: python app.py; open http://127.0.0.1:8000."""
import argparse
import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse, parse_qs
from bbb.pipeline import analyze, batch_analyze, uploaded_records, results_csv
from bbb.structures import index_status

ROOT = Path(__file__).parent
STATIC = ROOT / "web"
EXAMPLES = [
    {"compound_id": "Alprazolam", "smiles": "CC1=NN=C2CN=C(C3=CC=CC=C3)C3=C(C=CC(Cl)=C3)N12", "pKa": 5.01,
     "note": "Public paper example.", "source": "https://doi.org/10.1021/acs.jmedchem.9b01220.s003"},
    {"compound_id": "Caffeine", "smiles": "Cn1c(=O)c2c(ncn2C)n(C)c1=O", "pKa": None,
     "note": "Structure example; supply pKa to score."},
    {"compound_id": "Ethanol", "smiles": "CCO", "pKa": None, "note": "C2H6O; distinct from dimethyl ether."},
    {"compound_id": "Dimethyl ether", "smiles": "COC", "pKa": None, "note": "C2H6O; distinct from ethanol."},
]


class Handler(BaseHTTPRequestHandler):
    def local_request(self):
        allowed = {f"127.0.0.1:{self.server.server_port}", f"localhost:{self.server.server_port}"}
        host = self.headers.get("Host", "")
        origin = self.headers.get("Origin")
        if host not in allowed or (origin and origin not in {f"http://{h}" for h in allowed}):
            self.respond(403, {"error": "Local requests only."})
            return False
        return True

    def respond(self, status, content, kind="application/json; charset=utf-8", filename=None):
        data = json.dumps(content, ensure_ascii=False, allow_nan=False).encode() if kind.startswith("application/json") else content
        if isinstance(data, str):
            data = data.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", kind)
        self.send_header("Content-Length", str(len(data)))
        if filename:
            self.send_header("Content-Disposition", f'attachment; filename="{filename}"')
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Content-Security-Policy", "default-src 'self'; img-src 'self' blob: data:; style-src 'self'; script-src 'self'; connect-src 'self'; object-src 'none'; frame-ancestors 'none'")
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        if not self.local_request():
            return
        path = urlparse(self.path).path
        files = {"/": ("index.html", "text/html; charset=utf-8"),
                 "/app.js": ("app.js", "text/javascript; charset=utf-8"),
                 "/style.css": ("style.css", "text/css; charset=utf-8"),
                 "/template.csv": ("template.csv", "text/csv; charset=utf-8")}
        if path == "/api/examples":
            self.respond(200, EXAMPLES)
        elif path == "/api/structures-status":
            self.respond(200, index_status())
        elif path == "/api/health":
            from rdkit import rdBase
            self.respond(200, {"status": "ok", "rdkit": rdBase.rdkitVersion})
        elif path in files:
            filename, kind = files[path]
            self.respond(200, (STATIC/filename).read_bytes(), kind)
        else:
            self.respond(404, {"error": "Page not found."})

    def do_POST(self):
        # Local browser origin and Host checks keep private inputs on this server.
        if not self.local_request():
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if not 0 < length <= 24*1024*1024:
                raise ValueError("Request must be between 1 byte and 24 MB.")
            path = urlparse(self.path).path
            kind = self.headers.get("Content-Type", "")
            body = self.rfile.read(length)
            if "application/x-www-form-urlencoded" in kind and path in ("/api/export", "/api/download"):
                body = parse_qs(body.decode("utf-8")).get("payload", [""])[0]
            elif "application/json" not in kind:
                raise ValueError("JSON required.")
            payload = json.loads(body, parse_constant=lambda x: (_ for _ in ()).throw(ValueError("Non-finite JSON numbers are not allowed.")))
            if not isinstance(payload, dict):
                raise ValueError("Request must be a JSON object.")
            if path == "/api/analyze":
                self.respond(200, analyze(payload))
            elif path == "/api/batch":
                self.respond(200, batch_analyze(uploaded_records(payload)))
            elif path == "/api/export":
                rows = payload.get("results")
                if not isinstance(rows, list) or len(rows) > 10000 or not all(isinstance(r, dict) for r in rows):
                    raise ValueError("Invalid export data.")
                self.respond(200, "\ufeff"+results_csv(rows), "text/csv; charset=utf-8", "bbb_results.csv")
            elif path == "/api/download":
                data = analyze(payload.get("input"))
                if not data["structure"]:
                    raise ValueError("A valid structure is required for export.")
                if payload.get("format") == "svg":
                    self.respond(200, data["structure"]["svg"], "image/svg+xml; charset=utf-8", "molecule.svg")
                elif payload.get("format") == "mol":
                    self.respond(200, data["structure"]["molblock"], "chemical/x-mdl-molfile", "molecule.mol")
                else:
                    raise ValueError("Choose SVG or MOL.")
            else:
                self.respond(404, {"error": "Endpoint not found."})
        except (ValueError, TypeError, OverflowError) as error:
            self.respond(400, {"error": str(error)})
        except Exception:
            self.respond(500, {"error": "Processing failed. Check the file format."})

    def log_message(self, format, *args):
        # Log requests only; molecular payloads are never logged.
        print(f"{self.address_string()} {format % args}", flush=True)


def main():
    parser = argparse.ArgumentParser(description="BBB Score calculator")
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()
    server = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    print(f"BBB Score running at http://127.0.0.1:{args.port}", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
