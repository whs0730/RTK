import concurrent.futures
import socket
import subprocess
import sys
import time


def serve(listener):
    connection, _ = listener.accept()
    with connection:
        time.sleep(0.2)
        connection.sendall(bytes([0, 1]))
        time.sleep(0.05)
        connection.sendall(bytes([127, 128, 255]))
    connection, _ = listener.accept()
    with connection:
        connection.sendall(bytes([42]))


with socket.socket() as listener:
    listener.bind(("127.0.0.1", 0))
    listener.listen(2)
    listener.settimeout(10)
    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
        server = pool.submit(serve, listener)
        result = subprocess.run(
            [sys.argv[1], str(listener.getsockname()[1])],
            capture_output=True, text=True, timeout=10,
        )
        server.result(timeout=10)
    print(result.stdout, end="")
    if result.returncode:
        raise RuntimeError(result.stderr)
