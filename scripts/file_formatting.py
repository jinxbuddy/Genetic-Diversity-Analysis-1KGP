
with open("12_pop_1kg.geno") as f:
    data = f.readlines()
result = []
result2 = []

x = open("12_pop_reformatted.geno", "a")

for i in data:
    #print(i)
    for j in range(0,len(i)):
        if i[j] == "\n":
            x.write(i[j])

            
        else:
            x.write(i[j] + '\t')

x.close()


