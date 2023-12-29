import pathlib
import shutil
import tempfile
import json
import os
import subprocess
import signal
from concurrent.futures import ThreadPoolExecutor, Future
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

def evaluate(tree: dict, pool: ThreadPoolExecutor, prefix: str = "", idx: int | str = "") -> dict:
    """Recursively look for list of strings and execute each string as a python program.

    Note: key `_meta` will be ignored.
    It returns a dict with the same structure as tree except leaves that are lists of strings.
    Leaves that are lists of strings are replaced by a tuple (correct, total) describing the
    execution results.
    """
    result = {}
    for k, v in tree.items():
        if k == '_meta':
            result[k] = v
            continue
        if isinstance(v, list):
            if all(isinstance(elem, str) for elem in v):
                result[k] = [pool.submit(partial(_test, idx=idx or k, prefix=prefix), elem) for elem in v]
        elif isinstance(v, dict):
            result[k] = evaluate(v, pool, prefix=prefix, idx=idx or k)
        else:
            result[k] = v
    return result


def expand_future(node):
    """Recursively await and expand Future objects"""
    if isinstance(node, Future):
        return node.result()
    elif isinstance(node, dict):
        result = {}
        for k, v in node.items():
            result[k] = expand_future(v)
        return result
    elif isinstance(node, list):
        return [expand_future(v) for v in node]
    else:
        return node


def print_tree(node, depth):
    """Recursively print the tree"""
    indent = "  " * depth
    if isinstance(node, dict):
        for k, v in node.items():
            print(indent + k + ":")
            print_tree(v, depth+1)
    else:
        print(indent + str(node))


def _test(code: str, idx: int | str, prefix: str = "") -> tuple:
    folder_path = f"{prefix}/{path_prefix}/secret_{idx}"
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
            name = path.name.replace("_codes.json", "")
            res = evaluate(codes, pool, prefix=name)
            results[name] = res
    results = expand_future(results)
    print_tree(results, depth=0)

    json.dump(results, open("results.json", "w"))

if __name__ == '__main__':
    main()
