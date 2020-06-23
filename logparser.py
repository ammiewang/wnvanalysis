"""
parses a text/log file containing ecotypes into list of ecotypes
"""
def parser(inp):
    file = open(inp)
    ecotypes = []

    for line in file:
        if "Ecotype" in line and "[" in line and "]" in line:
            start = False
            ecotype = []
            sequence = ""
            for char in line:
                if char == '[':
                    start = True
                if start == True:
                    if char != '[' and char != ',' and char != ' ' and char != ']':
                        sequence += char
                    elif char == ' ':
                        ecotype.append(sequence)
                        sequence = ""
                    elif char == ']':
                        ecotype.append(sequence)
                        ecotypes.append(ecotype)
                        sequence = ""
                        break
    return ecotypes
