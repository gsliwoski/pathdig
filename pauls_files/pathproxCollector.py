
# Open ddg scores csv file
# search for pathprox score for each entry
# if found, append score and print to console
# done.

import argparse

parser = argparse.ArgumentParser(description='Collect pathprox scores and append to ddG scores')
parser.add_argument('ddGs', help="Inputfile containing all ddG scores from a ddGcollector.py run.")
parser.add_argument('pathprox', help="Folderstructure from a pathprox run prepared by pathproxPrepper.py")
parser.add_argument('--output', '-o', type=str, default='scoredVariants.csv', help='File to write out results.')
parser.add_argument('--report', '-r', type=str, default='pathproxFailureReport.txt', help='File to write a summary report of data collection.')

args = parser.parse_args()

ddGs = open(args.ddGs, 'r')

scores = open(args.output, 'w')
failureReport = open(args.report, 'w')

unpDict = {}
totalScores = 0
totalFails = 0

for line in ddGs:
    line = line[:-1]
    if 'structid' in line:
        scores.write(line + ",pathprox\n")
    else:
        linelist = line.split(',')
        structid = linelist[0]
        if structid not in unpDict:
            unpDict[structid] = [0, 0, []]
        chain = linelist[1]
        variant = linelist[10] + linelist[9] + linelist[11]
        hasScore = True
        try:
            summary = open(args.pathprox + structid + "/" + variant + "_out/" + structid + "_" + chain + "_D_summary.csv", 'r')
            for sumline in summary:
                if 'Kz_path' in sumline:
                    pass
                else:
                    pathproxScore = sumline.split("\t")[10]
                    if pathproxScore == 'nan':
                        raise Exception
            summary.close()
        except:
            hasScore = False
            unpDict[structid][1] += 1
            totalFails += 1
            unpDict[structid][2].append(variant)
        if hasScore:
            scores.write(line + "," + pathproxScore + "\n")
            unpDict[structid][0] += 1
            totalScores += 1
    
ddGs.close()
scores.close()

failureReport.write("Name\tScores\tFails")

totalStructs = 0
totalStructsWithFail = 0
totalStructsWithoutFails = 0
totalTotalFails = 0

for struct in unpDict:
    totalStructs += 1
    failureReport.write("\n" + struct + "\t" + str(unpDict[struct][0]) + "\t" + str(unpDict[struct][1]))
    if unpDict[struct][0] == 0:
        totalTotalFails += 1
    else:
        if unpDict[struct][1] == 0:
            totalStructsWithoutFails += 1
        else:
            totalStructsWithFail += 1
    
failureReport.write("\n" + "--------------------")
failureReport.write("\n" + "Sum:\t" + str(totalScores) + "\t" + str(totalFails) + "\n\n")
failureReport.write("\n" + "Total struct count: " + str(totalStructs))
failureReport.write("\n" + "Structs without fails: " + str(totalStructsWithoutFails))
failureReport.write("\n" + "Structs with both: " + str(totalStructsWithFail))
failureReport.write("\n" + "Failed structs: " + str(totalTotalFails) + "\n\n\n")

failureReport.write("Failed variants per structure:\n")
failureReport.write("[Mix ]: There are more variants for this structure, but they didn't fail.\n")
failureReport.write("[Fail]: These are all the variants of this structure, so all failed.\n")
for struct in unpDict:
    flag = ""
    if unpDict[struct][0] == 0:
        flag = "[FAIL]"
    else:
        if unpDict[struct][1] == 0:
            flag = ""
        else:
            flag = "[MIX ]"
    if flag != "":
        failureReport.write(flag + struct + ": " + str(unpDict[struct][2]) + "\n")

failureReport.close()
    