import pandas as pd

def hamming_distance(seq1, seq2): # calculate mismatches 
    if len(seq1) != len(seq2):
        raise ValueError("For Hamming distance, the length of sequence 1 must be equal to sequence 2!")
    else:
        for c1 in seq1:
            if c1.upper() not in ["A", "G", "T", "C"]:
                print(f"{seq1} is Invalid")
                return None
        for c2 in seq2:
            if c2.upper() not in ["A", "G", "T", "C"]:
                print(f"{seq2} is Invalid (contains {c2})")
                return None
        distance = 0
        for char1s, char2s in zip(seq1, seq2):
            if char1s != char2s:
                distance += 1
            
        # chars_1 = [char1 for char1 in char1s]
        # chars_2 = [char2 for char2 in char2s]
        return distance 

def l_matrix(seq1:str, seq2:str):
    """
    Aligns the sequence in a 2d matrix. 
    Diagonal is [i-1, j-1]
    if got [i, j-1], insertion -> j-1 means moving up 
    if got [i-1, j], deletion -> i-1 means moving left
    if got [i-1, j-1], match or substitution -> compare seq1[i-1] & seq2[j-1]
    --------------------------------
       | "" | A | G | T | C (seq2)
    --------------------------------
    "" |  0   1   2   3   4        -> at A,I'm trying to align A with an empty string->needs deletion(moves left)& 0+1 
    --------------------------------
    A  |  1   0   1   2            -> at A, im trying to align A with an empty string->needs insertion(moves up)& 0+1 
    --------------------------------
    G  |  2   

    Examples:
    for (A,G): diagonal mismatch: 1+1 = 2;
    Top: 2+1 = 3
    Left: 0+1 = 1
    min(2,3,1) = 1 ->  insertion
    for (A,T):
    diagonal: 2+1 = 3
    Top: 3+1 = 4
    Left: 1+1 = 2
    min(3,4,2) = 2 -> insertion
    conclusion: like trying to align AGTC with AG, i.e. TC has been inserted
    good visualization: https://www.youtube.com/watch?v=lBU2rZe-whk  
    """
    # creating a blank matrix 
    rows = len(seq1)+1
    cols = len(seq2)+1
    matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    # matrix = [[0] * cols] * rows # don't use this
    # matrix[0][1] = 1 # modifying this would change column 1 across every row 
    # i.e.[[0, 1, 0, 0, 0], [0, 1, 0, 0, 0], [0, 1, 0, 0, 0], [0, 1, 0, 0, 0], [0, 1, 0, 0, 0]]
    for i in range(rows):
        matrix[i][0] = i
    for j in range(cols):
        matrix[0][j] = j 
    # for row in matrix:
    #     print(row) # black canvas: [[0, 1, 2, 3, 4], [1, 0, 0, 0, 0], [2, 0, 0, 0, 0], [3, 0, 0, 0, 0], [4, 0, 0, 0, 0]]
        # 0's are filled with the alignment 
    
    # the nested loop 
    for i in range(1, rows):
        for j in range(1, cols):
            cost = 0 if seq1[i-1] == seq2[j-1] else 1 
            deletion = matrix[i-1][j] + 1 # +1 for penalty for insertion + deletion 
            substitution = matrix[i-1][j-1] + cost
            insertion = matrix[i][j-1] +1 
            matrix[i][j] = min(deletion, substitution, insertion)
            # if matrix[i][j] == deletion:
            #     print(f"deletion")

    # traceback
    i = rows-1 # moving left
    j = cols-1 # moving up
    operations = []
    top_row = []
    match_row = []
    bot_row = []
    while i > 0 or j > 0:
        # must be defined here so that its initialised for other conditions 
        current = matrix[i][j] # check which neighbour produces the current cell 
        # top_row.append(seq1[i-1]) # cannot put them here, because when j stays the same (deletion), this is still appended on the 2nd loop,
        # causing characters to be duplicated 
        
        
        # bot_row.append(seq2[j-1])
        if i > 0 and j > 0: 
            cost = 0 if seq1[i-1] == seq2[j-1] else 1 # must be recalculated otherwise conatins the very last cell of nested loop above 
            if current == matrix[i-1][j-1] + cost: # substitution condition
                top_row.append(seq1[i-1])
                bot_row.append(seq2[j-1])
                if seq1[i-1] == seq2[j-1]:
                    operations.append(f"Match of {seq1[i-1]} at position {i-1}")
                    match_row.append("|")
                else:
                    operations.append(f"Substitution of {seq1[i-1]} in sequence 1 to {seq2[j-1]} in sequence 2 at position {i-1}")
                    match_row.append("*")
                # has to be defined here so that it only moves diagonally if the path is diagonal 
                i -= 1
                j -= 1 
                continue # must be here so only runs if the substition condition is True 
        if i > 0: # if deletion j wouldve arrived at 0 faster than i 
            if current == matrix[i-1][j] + 1:
                operations.append(f"Deletion of {seq1[i-1]} in sequence 1 at position {i-1}")
                top_row.append(seq1[i-1])
                bot_row.append("-")
                match_row.append(" ")
                i -= 1
                continue 
        if j > 0:
            if current == matrix[i][j-1] + 1:
                operations.append(f"Insertion of {seq2[j-1]} in sequence 2 at position {j-1}")
                top_row.append("-")
                match_row.append(" ")
                bot_row.append(seq2[j-1])
                j -= 1
                continue 

    operations.reverse() # reverse to show operations from start to end cz currently it is from end to start (traceback) 
    bot_row.reverse()
    top_row.reverse()
    match_row.reverse()
    return matrix[rows-1][cols-1], operations, top_row, match_row, bot_row, 

def main():
    seq_1 = str(input("Sequence 1: "))
    seq_2 = str(input("Sequence 2: "))
    print(f"Matching {seq_2} with {seq_1}")
    try:
        dist = hamming_distance(seq_1, seq_2)
        print(f"hamming_distance: {dist}")
    except ValueError as err:
        print(f"Note: {err}")

    
    matrix, operations, top, match, bot = l_matrix(seq_1, seq_2)
    print(f"levenshtein distance: {matrix}")
    print(f"Sequence 1: {"".join(top)}")
    print(f"            {"".join(match)}")
    print(f"Sequence 2: {"".join(bot)}")
    print(f"- means insertion or deletion; * means substitution")
    print(f"Operations: {operations}")
    
    #print(f"return: {matrix}")

    #print(f"{matrix[0][1]}")

if __name__== "__main__":
    main()