import re
import copy
from configuration import indel, similarity_matrix


def init():
    global sm
    sm = {}

    for yi, y in enumerate("AGCT"):
        sm[y] = {}
        for xi, x in enumerate("AGCT"):
            sm[y][x] = similarity_matrix[yi][xi]


def main():
    init()

    dna1 = input("Sequence1< ").upper()
    dna2 = input("Sequence2< ").upper()
    dna3 = input("Sequence3< ").upper()
    dna1 = re.sub(' ', '', dna1)
    dna2 = re.sub(' ', '', dna2)
    dna3 = re.sub(' ', '', dna3)

    table = create_table(dna1, dna2, dna3)


def create_table(dna1, dna2, dna3):
    # create table
    table = []
    for y in range(len(dna2) + 1):
        table.append([])
        for z in range(len(dna3) + 1):
            table[y].append([])
            for x in range(len(dna1) + 1):
                table[y][z].append({"value": 0, "arrows": [False, False, False, False]})
            #               ^ A cell of the Table.   Left   Diag.  Up

    # fill first row
    for i in range(len(dna1) + 1):
        table[0][0][i]["value"] = 5 * i
        table[0][0][i]["arrows"][0] = True

    # fill first column
    for i in range(len(dna2) + 1):
        table[i][0][0]["value"] = 5 * i
        table[i][0][0]["arrows"][2] = True

    for i in range(len(dna3) + 1):
        table[0][i][0]["value"] = 5 * i
        table[0][i][0]["arrows"][2] = True

    # while filling first row and column, we allowed this mistake for the sake
    # of simplicity but we need to fix that
    table[0][0][0]["arrows"] = [False, False, False, False]

    # calculate the rest
    Cost = 0
    for i2, b2 in enumerate(dna2, start=0):
        for i3, b3 in enumerate(dna3, start=0):
            for i1, b1 in enumerate(dna1, start=0):
                score1 = table[i2][i3 - 1][i1 - 1]["value"] + indel
                score2 = table[i2 - 1][i3 - 1][i1 - 1]["value"] + calculate_scoreNew(dna1[i1], dna2[i2], dna3[i3])
                score3 = table[i2 - 1][i3 - 1][i1]["value"] + indel
                score4 = table[i2 - 1][i3][i1 - 1]["value"] + indel
                min_score = min(score1, score2, score3, score4)
                table[i2][i3][i1]["value"] = min_score
                Cost = min_score

                if score1 == min_score:
                    table[i2][i3][i1]["arrows"][0] = True
                if score2 == min_score:
                    table[i2][i3][i1]["arrows"][1] = True
                if score3 == min_score:
                    table[i2][i3][i1]["arrows"][2] = True
                if score4 == min_score:
                    table[i2][i3][i1]["arrows"][3] = True
    print("\t Optimal Cost: %s\t" % Cost)
    return table


def print_table(dna1, dna2, table):
    # Numbers Table
    print("Scores:")
    print("\t      ", end="")
    for i in dna1:
        print(i, end="  ")
    print()

    for iy, y in enumerate(table):
        print("\t%s" % (" " + dna2)[iy], end=" ")
        for x in y:
            print("% d" % x["value"], end=" ")
        print()

    print()

    # Arrows Table
    print("Arrows: (value = 1 if left + 2 if diagonal + 4 if up)")
    print("\t      ", end="")
    for i in dna1:
        print(i, end="  ")
    print()

    for iy, y in enumerate(table):
        print("\t%s" % (" " + dna2)[iy], end=" ")
        for x in y:
            summe = x["arrows"][0] + x["arrows"][1] * 2 + x["arrows"][2] * 4
            print("% d" % summe, end=" ")
        print()

    print("\n")


def find_alignments(dna1, dna2, dna3, table):
    alignments = []

    # trace_r: trace arrows back to their origin
    # !recursive function!
    def trace_r(y, x, z, alignment):
        if table[y - 1][x - 1][z - 1]["arrows"][0]:  # LEFT
            alignment_c = copy.deepcopy(alignment)  # create a branch
            alignment_c["DNAs"][0] += dna1[x - 1]
            alignment_c["DNAs"][1] += "-"
            alignment_c["DNAs"][2] += dna3[z - 1]
            alignment_c["score"] += table[y - 1][x - 1]["value"]
            trace_r(y, x - 1, z - 1, alignment_c)
            del alignment_c

        if table[y - 1][x - 1][z - 1]["arrows"][1]:  # DIAGONAL
            alignment_c = copy.deepcopy(alignment)  # create a branch
            alignment_c["DNAs"][0] += dna1[x - 1]
            alignment_c["DNAs"][1] += dna2[y - 1]
            alignment_c["DNAs"][2] += dna3[z - 1]
            alignment_c["score"] += table[y - 1][x - 1][z - 1]["value"]
            trace_r(y - 1, x - 1, z - 1, alignment_c)
            del alignment_c

        if table[y - 1][x - 1][z - 1]["arrows"][2]:  # UP
            alignment_c = copy.deepcopy(alignment)  # create a branch
            alignment_c["DNAs"][0] += "-"
            alignment_c["DNAs"][1] += dna2[y - 1]
            alignment_c["DNAs"][2] += dna3[z - 1]
            alignment_c["score"] += table[y - 1][x - 1]["value"]
            trace_r(y - 1, x, z - 1, alignment_c)
            del alignment_c

        # origin :)
        if y - 1 == -len(table) and x - 1 == -len(table[0]) and z - 1 == -len(table[1]):
            alignment["DNAs"] = alignment["DNAs"][0][::-1], alignment["DNAs"][1][::-1], alignment["DNAs"][2][::-1]
            score = calculate_score(alignment["DNAs"][0], alignment["DNAs"][1], alignment["DNAs"][2],
                                    alignment["score"])
            alignments.append({"DNAs": alignment["DNAs"], "score": score})
        trace_r(0, 0, 0, {"DNAs": ["", "", ""], "score": 0})

    return alignments


def print_alignments(alignments):
    print("%d Alignments:" % (len(alignments)))
    for a in alignments:
        print("\t%s\tScore:" % (a["DNAs"][0]))
        print("\t%s\t%d" % (a["DNAs"][1], a["score"]))
        print()


def calculate_score(alg1, alg2, score):
    score = 0
    for i in range(len(alg1)):
        if alg1[i] == '-' or alg2[i] == '-':
            score += indel
        elif alg1[i] == alg2[i]:
            score += 0
        elif alg1[i] == 'T' and (alg2[i] == 'A' or alg2[i] == 'G'):
            score += 5
        elif alg1[i] == 'A' and (alg2[i] == 'T' or alg2[i] == 'C'):
            score += 5
        elif alg1[i] == 'C' and (alg2[i] == 'A' or alg2[i] == 'G'):
            score += 5
        elif alg1[i] == 'G' and (alg2[i] == 'C' or alg2[i] == 'T'):
            score += 5
        elif alg1[i] == ' ' and alg2 == ' ':
            score = score
        else:
            score += 2

    return score


def calculate_scoreNew(alg1, alg2, alg3):
    score = 0

    if (alg1 == alg2) and (alg3 == alg1):
        score += 0

    elif alg1 == 'T' and (alg2 == 'A' or alg2 == 'G') and (alg3 == 'A' or alg3 == 'G'):
        score += 5

    elif alg1 == 'A' and (alg2 == 'T' or alg2 == 'C') and (alg3 == 'T' or alg3 == 'C'):
        score += 5

    elif alg1 == 'C' and (alg2 == 'A' or alg2 == 'G') and (alg3 == 'A' or alg3 == 'G'):
        score += 5

    elif alg1 == 'G' and (alg2 == 'C' or alg2 == 'T') and (alg3 == 'C' or alg3 == 'T'):
        score += 5

    elif alg1 == ' ' and alg2 == ' ':
        score = score

    else:
        score += 2

    return score


if __name__ == "__main__":
    main()
