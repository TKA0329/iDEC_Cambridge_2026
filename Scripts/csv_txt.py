import pandas as pd

df = pd.read_csv('sequences.csv')
with open('output.txt', 'w') as f:
    for i, seq in enumerate(df['sequence'], 1):
        f.write(f'>{i}\n{seq}\n')