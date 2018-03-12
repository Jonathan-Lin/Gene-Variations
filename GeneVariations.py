import csv
import random
import statistics
import math
import matplotlib.pyplot as plt

#Generate variants based on frequency and generate plot from that

#using a probability, if random number is less than, true, if greater than, false
def decision(probability):
    probable = float(probability);
    randNum = random.uniform(0, 1);
    
    if randNum < probable:
        return 1
    else:
        randNum = random.uniform(0, 1);
        if randNum < probable:
            return 1
        else:
            return 0

#calculate the decisions based on list of frequencies passed 
def calculator(list_of_freq):
    variants = 0
    for i in range(len(list_of_freq)):
        if(decision(list_of_freq[i]) == 1):
            variants += 1
    return variants

#Used to label the number above each bar in the chart
def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        plt.text(rect.get_x()+rect.get_width()/2., 1.0*height, '%d'%int(height),
                ha='center', va='bottom')

#plot the list of variants resulted from calculations and generate histogram
def plot(list_of_results, possible, csv_file, variant_type ):
    
    plot_data = [0] * possible
    for x in list_of_results:
        plot_data[x] = plot_data[x] + 1

    rect = plt.bar(range(len(plot_data)), plot_data, align='center', alpha=0.5)

    plt.ylabel('Number of Occurances with 10000 Runs')
    plt.xlabel('Number of Variants Seen')
    if (variant_type == "1"):
        plt.title('Distribution of # of Variants Seen for ' + str(csv_file) + ' Gene')
    elif (variant_type == "2"):
        plt.title('Distribution of # of Rare Variants Seen for ' + str(csv_file) + ' Gene')
        
    plt.xticks(range(0, possible, 1))
    
    autolabel(rect)
    
    plt.show()

#Calculate average variants, stdev, for rare variants and variants
def average_variants(start, end, race, csv_file):
    with open(csv_file,'rt') as file:
        reader = csv.reader(file)
        my_list = []
        firstLine = True

        #For each row, add the information we will be using, sorry for hard coded numbers here
        for row in reader:
            if firstLine:
                print("You are looking at", row[3 * int(race) + 12 ],row[3 * int(race) + 13 ])
                firstLine = False
                continue
            #This calculates which column for each race
            position = (3 * int(race)) + 11 
            my_list.append([row[1], row[9], row[position], row[position+1], row[position+2], row[position+3]])

    in_range = False
    frequency_list = []
    num_introns = 0

    #generate frequencies from allele number and variants counted
    for row in my_list:
        #start adding when in range, or if -7 then means look at entire gene
        if (row[0] == start or start==-7):
            in_range = True

        #This part will skip over intronic and exonic parts we will not look at
        if(row[1]=="5\' UTR" or row[1]=="intron" or row[1]=="3\' UTR"):
            num_introns += 1
            continue
        if(row[1]=="upstream gene" or row[1]=="downstream gene" or row[1]=="non coding transcript exon"):
            num_introns += 1
            continue
        
        #add frequency if in range
        if(in_range):
            try:
                if(race==0):
                    frequency_list.append((float(row[2])/float(row[3])))
                else:
                    frequency_list.append((float(row[3])/float(row[4])))
            except ZeroDivisionError:
                print("Had 0 cases of variant", row[0],row[1],row[2], row[3], row[4])
                
        #stop when out of range
        if row[0] == end:
            in_range = False
                
    print("Frequency List")
    print(len(frequency_list))
    
    rare_frequency_list = []

    #THIS IS WHERE TO CHANGE WHAT FREQUENCY CUT OFF TO LOOK AT CURRENTLY < 1%
    for frequency in frequency_list:
        if ( frequency < .01 ):
            rare_frequency_list.append(frequency)

    print("\nRare Frequency List")
    print(len(rare_frequency_list))
    print()
            
    total = 0;
    rare_total = 0
    runs = 10000;
    var_each_run_list = []
    rare_var_each_run = []
    #Run 10,000 times and keep track of all runs here
    for x in range(runs):
        runCount = calculator(frequency_list)
        var_each_run_list.append(runCount)
        rareRunCount = calculator(rare_frequency_list)
        rare_var_each_run.append(rareRunCount)

        total += runCount
        rare_total += rareRunCount
    
    average = (float(total) / float(runs))
    rare_average = (float(rare_total) / float(runs))
    print( "Total Variants Over 10000 Runs: ", total )
    print( "Total Rare Variants over 10000 Runs ", rare_total)
    print( "Average of Variants over 10000 Runs: " + str(average))
    print( "Average Rare Variants over 10000 Runs: " + str(rare_average))
    print( "Standard Deviation of Population: " + str(statistics.pstdev(var_each_run_list)))
    print( "Standard Deviation of Rare Population: " + str(statistics.pstdev(rare_var_each_run)))

    #Calculate Confidence Interval
    #sqrt_sample_size = math.sqrt(len(frequency_list))
    #first_step = (statistics.pstdev(var_each_run_list) / sqrt_sample_size)
    #deviation = (1.96 * first_step)
    #last_step_pos = (average + deviation)
    #last_step_neg = (average - deviation)
    #print ("Confidence Interval of 95%: Low: " + str(last_step_neg) + " High: " + str(last_step_pos))
    #print()

    #Ask if you want to plot
    choice = input("Plot 1. All variants 2. Rare Variants 3. None\n")
    if( choice == "1" ):
        plot(var_each_run_list, len(frequency_list), csv_file, '1')
    elif( choice == "2"):
        plot(rare_var_each_run, len(frequency_list), csv_file, '2')

        
#Choices for gene within aHUS Panel
def ahus_choice():
    gene_choice = input("What Gene?\n1.CFH\n2.CFHR1\n3.CFHR3\n4.CFHR4\n5.CFHR5\n6.C3\n7.CFB\n8.CFI\n9.THBD\n10.DGRE\n11.CD46\n12.PLG\n")

    if gene_choice == "1":
        gene_choice = 'CFH.csv'
    elif gene_choice == "2":
        gene_choice = 'CFHR1.csv'
    elif gene_choice == "3":
        gene_choice = 'CFHR3.csv'
    elif gene_choice == "4":
        gene_choice = 'CFHR4.csv'
    elif gene_choice == "5":
        gene_choice = 'CFHR5.csv'
    elif gene_choice == "6":
        gene_choice = 'C3.csv'
    elif gene_choice == "7":
        gene_choice = 'CFB.csv'
    elif gene_choice == "8":
        gene_choice = 'CFI.csv'
    elif gene_choice == "9":
        gene_choice = 'THBD.csv'
    elif gene_choice == "10":
        gene_choice = 'DGKE.csv'
    elif gene_choice == "11":
        gene_choice = 'CD46.csv'
    elif gene_choice == "12":
        gene_choice = 'PLG.csv'
    else:
        gene_choice = "Invalid"
    return gene_choice

#The Script
#Ask which gene or panel you want
panel_choice = input("Which Panel? \n1.aHus\n2.MTHFR\n3.VWF\n4.ADAMTS13\n5.COL4A3\n6.COL4A4\n7.COL4A5\n")

#Gene choice reflects which number is inputted
#ADD NEW GENE NAMES HERE
if panel_choice == "1": gene_choice = ahus_choice()
elif panel_choice == "2": gene_choice = "MTHFR.csv"
elif panel_choice == "3": gene_choice = "VWF.csv"
elif panel_choice == "4": gene_choice = "ADAMTS13.csv"
elif panel_choice == "5": gene_choice = "COL4A3.csv"
elif panel_choice == "6": gene_choice = "COL4A4.csv"
elif panel_choice == "7": gene_choice = "COL4A5.csv"
else: gene_choice = "Invalid"
    
if gene_choice != "Invalid":
    choice = input("What do you want to run? \n1. By Race 2. By Range\n")

    if (choice == '1'):
        print ("You entered", choice)
        raceChoice = input("What Race? African, East Asian, European(Non-Finnish), Finnish, Latino, Other, South Asian\n")       
    elif ( choice == '2'):
        raceChoice = 0
    else: print("Not a Valid Gene Choice")

    if (choice == '1'):
            range_choice = input("Do you want to specify a range? y or n\n")
            if (range_choice=='n'):
                startInput = -7
                endInput = -7
            elif(range_choice=='y'):
                startInput = input("Starting Position?\n")
                endInput = input("Ending Position?\n")
    else:
        startInput = input("Starting Position?\n")
        endInput = input("Ending Position?\n")           
    
    average_variants(startInput, endInput, raceChoice, gene_choice)
    
else: print("Invalid Choice")
 

    
