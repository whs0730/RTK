import csv
from pathlib import Path
import signal
import socket
import subprocess
import sys
import tempfile
import time

binary = str(Path(sys.argv[1]).resolve())


def run_case(stop_signal=None):
    with tempfile.TemporaryDirectory(prefix="rtk-shutdown-", dir=Path.cwd()) as temporary:
        directory = Path(temporary)
        with socket.socket() as base, socket.socket() as rover:
            for listener in (base, rover):
                listener.bind(("127.0.0.1", 0))
                listener.listen(1)
                listener.settimeout(10)
            with (directory / "run.log").open("w") as log:
                process = subprocess.Popen(
                    [binary, "--stream-rooftop",
                     "127.0.0.1", str(base.getsockname()[1]),
                     "127.0.0.1", str(rover.getsockname()[1])],
                    cwd=directory, stdout=log, stderr=subprocess.STDOUT,
                )
                try:
                    base_peer, _ = base.accept()
                    rover_peer, _ = rover.accept()
                    with base_peer, rover_peer:
                        deadline = time.monotonic() + 10
                        while "Realtime processing started." not in (directory / "run.log").read_text():
                            if process.poll() is not None or time.monotonic() > deadline:
                                raise RuntimeError("stream did not start")
                            time.sleep(0.02)
                        if stop_signal is None:
                            base_peer.shutdown(socket.SHUT_WR)
                        else:
                            process.send_signal(stop_signal)
                        result = process.wait(timeout=5)
                finally:
                    if process.poll() is None:
                        process.kill()
                        process.wait()
            text = (directory / "run.log").read_text()
            if stop_signal is None:
                assert result != 0 and "realtime input failed" in text, text
            else:
                assert result == 0 and "Finished." in text, text
                assert f"Stop signal: {stop_signal.value}" in text, text
            assert "Synchronized epochs: 0" in text, text
            for stem in ("rtk_float", "rtk_solution"):
                with (directory / f"{stem}_realtime_rooftop.csv").open() as output:
                    rows = list(csv.reader(output))
                assert len(rows) == 1 and len(rows[0]) > 10, rows
            assert len(list(directory.glob("*.txt"))) == 7
            print(f"{stop_signal or 'EOF'}: exit={result}, CSV headers and final statistics saved")


for stop_signal in (signal.SIGINT, signal.SIGTERM, None):
    run_case(stop_signal)
