#s30995
#05.05.2026
#Random DNA sequence generator in FASTA format

import random


def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000) -> int:
    while True:
        value = input(prompt)
        try:
            num = int(value)
            if min_val <= num <= max_val:
                return num
        except ValueError:
            pass
        print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")


def get_valid_id() -> str:
    while True:
        seq_id = input("Enter sequence ID: ")
        if seq_id.strip() == "" or any(c.isspace() for c in seq_id):
            print("Error: ID must not contain whitespace.")
        else:
            return seq_id


def get_distribution() -> dict:
    print("Enter nucleotide percentages (must sum to 100):")
    while True:
        try:
            a = float(input("A: "))
            c = float(input("C: "))
            g = float(input("G: "))
            t = float(input("T: "))
            total = a + c + g + t
            if abs(total - 100.0) < 1e-6:
                return {"A": a, "C": c, "G": g, "T": t}
        except ValueError:
            pass
        print("Error: values must be numbers summing to 100.")


def generate_sequence(length: int, distribution: dict = None) -> str:
    nucleotides = ["A", "C", "G", "T"]

    if distribution:
        weights = [distribution[n] for n in nucleotides]
        return "".join(random.choices(nucleotides, weights=weights, k=length))
    else:
        return "".join(random.choice(nucleotides) for _ in range(length))


def calculate_stats(sequence: str) -> dict:
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    total = 0

    for char in sequence:
        if char in counts:
            counts[char] += 1
            total += 1

    stats = {}
    for n in counts:
        stats[n] = (counts[n] / total * 100) if total > 0 else 0.0

    gc_ratio = ((counts["G"] + counts["C"]) / total * 100) if total > 0 else 0.0
    stats["gc_ratio_A"] = gc_ratio

    return stats


def insert_name(sequence: str, name: str) -> str:
    pos = random.randint(0, len(sequence))
    return sequence[:pos] + name.lower() + sequence[pos:]


def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:
    header = f">{seq_id}"
    if description:
        header += f" {description}"

    lines = [sequence[i:i + line_width] for i in range(0, len(sequence), line_width)]
    return header + "\n" + "\n".join(lines) + "\n# EOF_1\n"


def find_motif(sequence: str, motif: str) -> list:
    positions = []
    for i in range(len(sequence) - len(motif) + 1):
        if sequence[i:i + len(motif)] == motif:
            positions.append(i + 1)
    return positions


def reverse_complement(sequence: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp.get(base, base) for base in reversed(sequence))


def transcribe(sequence: str) -> str:
    return sequence.replace("T", "U")


def main():
    length = validate_positive_int("Enter sequence length: ")
    seq_id = get_valid_id()
    description = input("Enter a description of the sequence: ")
    name = input("Enter your name: ")

    use_dist = input("Custom nucleotide distribution? (y/n): ").lower()
    distribution = get_distribution() if use_dist == "y" else None

    sequence = generate_sequence(length, distribution)

    stats = calculate_stats(sequence)

    sequence_with_name = insert_name(sequence, name)

    fasta_content = format_fasta(seq_id, description, sequence_with_name)

    filename = f"{seq_id}.fasta"
    with open(filename, "w") as f:
        f.write(fasta_content)

    print(f"\nSequence saved to file: {filename}\n")

    print(f"Sequence statistics (n={length}):")
    for n in ["A", "C", "G", "T"]:
        print(f"{n}: {stats[n]:.2f}%")
    print(f"GC-content: {stats['gc_ratio_A']:.2f}%")

    motif = input("\nEnter motif to search (or empty to skip): ").upper()
    if motif:
        positions = find_motif(sequence, motif)
        print("Motif positions:", positions if positions else "Not found")

    rev_comp = reverse_complement(sequence)
    print("\nReverse complement (first 60 nt):")
    print(rev_comp[:60])

    mrna = transcribe(sequence)
    print("\nmRNA (first 60 nt):")
    print(mrna[:60])


if __name__ == "__main__":
    main()