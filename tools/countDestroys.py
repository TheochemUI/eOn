import re
from pathlib import Path

# This helps find leaks, since the potential calls will not be correct if there's an unpaired new/delete..

def extract_sum_calls(log_file_path, regex, count_lines_only=False):
    pattern = re.compile(regex)
    log_content = Path(log_file_path).read_text()
    matches = pattern.findall(log_content)
    if count_lines_only:
        total_calls = len(matches)
    else:
        total_calls = sum(map(int, matches))
    return total_calls


log_file_path = "_potcalls.log"
destroyed_calls = extract_sum_calls(log_file_path, r"destroyed after (\d+) calls")
all_calls = extract_sum_calls(log_file_path, r"so far", True)
print(f"Destructor counts: {destroyed_calls}")
print(f"All counts: {all_calls}")
