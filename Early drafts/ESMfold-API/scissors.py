with open('protein.fasta', 'r') as f:
    content = f.read().splitlines()
    sequence = ''.join(content[1:])

positions = input('Enter positions separated by commas: ')
positions = [int(pos) for pos in positions.split(',')]

with open('output.txt', 'w') as f:
    for pos in positions:
        start = max(0, pos - 200)
        end = min(len(sequence), pos + 199)
        f.write(sequence[start:end] + '\n')