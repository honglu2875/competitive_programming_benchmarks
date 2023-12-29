import requests
import pathlib
import shutil
import tempfile
import json
import os
import subprocess
import signal
from concurrent.futures import ThreadPoolExecutor
from functools import partial


path_prefix = "dataset"

def compare_files_ignore_empty_lines(file1, file2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        lines1 = [line.strip() for line in f1.readlines() if line.strip()]
        lines2 = [line.strip() for line in f2.readlines() if line.strip()]

    return lines1 == lines2

def run_with_timeout(cmd, timeout_sec, cwd=None):
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, cwd=cwd)
    try:
        output, error = process.communicate(timeout=timeout_sec)
        return process.returncode, output.decode('utf-8'), error.decode('utf-8')
    except subprocess.TimeoutExpired:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        return -1, None, None

def _test(code: str, idx: int | str):
    folder_path = f"{path_prefix}/secret_{idx}"
    # Loop through all filenames under the folder
    count = 0
    total = 0
    timeout_sec = 1
    for filename in os.listdir(folder_path):
        with tempfile.TemporaryDirectory() as temp_dir:
            if filename.endswith('.in'):
                total += 1
                input_tempfile_path = os.path.join(temp_dir, 'input.txt')
                output_tempfile_path = os.path.join(temp_dir, 'output.txt')
                program_tempfile_path = os.path.join(temp_dir, 'temp_program.py')

                input_file = os.path.join(folder_path, filename)
                output_file = os.path.join(folder_path, filename.replace('.in', '.ans'))

                # Copy the content of <filename>.in to input temporary file
                if os.path.exists(input_tempfile_path):
                    os.remove(input_tempfile_path)
                shutil.copy(input_file, input_tempfile_path)

                # Run the Python program to create output.txt
                python_code = code

                # Write the Python code to a temporary .py file
                with open(program_tempfile_path, 'w') as program_file:
                    program_file.write(python_code)

                # Run the temporary Python program to create output.txt with a timeout
                cmd = ['python', program_tempfile_path]
                return_code, _, _ = run_with_timeout(cmd, timeout_sec, cwd=temp_dir)

                # Compare the content of output.txt and <filename>.ans
                if return_code == 0:
                    if os.path.exists(output_tempfile_path):
                        if compare_files_ignore_empty_lines(output_tempfile_path, output_file):
                            count += 1
    return count, total


def main():
    collected = list(pathlib.Path(".").glob("*_codes.json"))
    print("Collected:")
    for i, path in enumerate(collected):
        print(f"  {i+1}. {path}")

    results = {}
    for path in collected:
        res = {}
        codes = json.load(open(str(path)))
        with ThreadPoolExecutor() as pool:
            for k, texts in codes.items():
                res[k] = []
                for text in texts:
                    res[k].append(pool.submit(partial(_test, idx=k), text))

            for k in codes:
                res[k] = [f.result() for f in results[k]]
                print(res[k])

        results[path.name.replace("_codes.json", "")] = res

    json.dump(results, open("results.json", "w"))

if __name__ == '__main__':
    main()
