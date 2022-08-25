import re
def ReadROH(inputt):
    with open(inputt, 'r') as T:
        for line in T:
            if line.startswith("## INFO:"):
                line = line.split(":")[1]
                line = float(line.split("Mb")[0].strip())
                return line
                # return float(re.findall(r'\d\d', line)[0])