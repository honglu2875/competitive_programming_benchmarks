import requests
import pathlib
import json
import re


path_prefix = "dataset"
instr = "### Instruction: Please write a python program that reads an input file input.txt and put the output in output.txt according to the following problem description. Make sure to wrap your python program inside ```python tag."
suffix = "### Response:"
pattern = re.compile(r'```python(.*?)```', re.DOTALL)


def extract_python_code(markdown_text):
    # Use regex to extract Python code within ```python``` tags
    idx = markdown_text.find("### Response:")
    markdown_text = markdown_text[idx+len(suffix):]
    match = pattern.search(markdown_text)
    if match:
        return match.group(1).strip()
    else:
        return None

def _get_order(s: str):
    match = re.search(r"problem_([\d]+)\.txt", s)
    assert match is not None
    return match.group(1)

def _contain_tests(p: pathlib.Path, order: str):
    return (p / f"secret_{order}").is_dir()

def _ensure_server():
    url = "http://localhost:8000/health"
    try:
        r = requests.get(url)
        if r.ok:
            return True
    except:
        pass
    return False

def _rename(name: str):
    return name.split(".")[0].replace("-", "_")


def main():
    collected = [p for p in pathlib.Path(".").iterdir() if (p / path_prefix).is_dir()]

    if not collected:
        print("No directories found that contains `dataset` folder.")
        exit(0)

    print("Collected:")
    for i, name in enumerate(collected):
        print(f"  {i+1}. {name}")

    for folder in collected:
        logs = {}
        problem_files = list((folder / path_prefix).glob("problem_*.txt"))
        problem_files = [p for p in problem_files if _contain_tests(p.parent, _get_order(p.name))]

        if not problem_files:
            print(f"No problems found in {folder}.")
            continue
        print(f"Found {len(problem_files)} problems in {folder}.")

        if not _ensure_server():
            print("Server is not running. Please run `python -m vllm.entrypoints.api_server --model ...` first.")
            continue

        orders = [_get_order(p.name) for p in problem_files]
        for i, f in zip(orders, problem_files):
            with open(f, "r") as fl:
                desc = fl.read()
            prompt = f"{instr}\n{desc}\n{suffix}"
            payload = {
                "prompt": prompt,
                "use_beam_search": False,
                "temperature": 0.2,
                "max_tokens": 800,
                "n": 10,
            }
            res = requests.post("http://localhost:8000/generate", json=payload)
            texts = res.json()["text"]
            logs[i] = [extract_python_code(s) for s in texts]

        open(f"{_rename(folder.name)}_codes.json", "w").write(json.dumps(logs))


if __name__ == "__main__":
    main()

