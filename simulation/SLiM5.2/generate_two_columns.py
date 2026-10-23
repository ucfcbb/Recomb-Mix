import sys


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    input_strings = sys.argv[3:]
    with open(input_file, "r") as f:
        lines = [line.strip() for line in f if line.strip() != ""]
    prefix_to_string = {}
    next_string_idx = 0
    prefix_counter = {}
    output_lines = []
    for line in lines:
        prefix, _, _ = line.partition(":")
        if prefix not in prefix_to_string:
            if next_string_idx >= len(input_strings):
                raise ValueError(f"Found prefix '{prefix}' but ran out of input strings (only {len(input_strings)} provided).")
            prefix_to_string[prefix] = input_strings[next_string_idx]
            next_string_idx += 1
            prefix_counter[prefix] = 0
        counter = prefix_counter[prefix]
        prefix_counter[prefix] += 1
        # a = counter // 2
        # b = counter % 2
        # label = f"{prefix_to_string[prefix]}_{a}_{b}"
        label = f"{prefix_to_string[prefix]}_{counter}"
        output_lines.append(f"{line}\t{label}")
    with open(output_file, "w") as f:
        f.write("\n".join(output_lines) + "\n")
    sys.exit(0)
