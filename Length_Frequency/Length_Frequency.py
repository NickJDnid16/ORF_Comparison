


Min_0 = open('/home/nick/Git/Open_Reading_Frame_Comparison/STORF/Pseudo/STORF_Pseudo.txt', mode='rb')
Lengths = []
m10 = []
m50 = []
m100 = []
m200 = []
m350 = []
m500 = []
m650 = []
m800 = []
m1000 = []

for line in Min_0:
    if "+" in line or "-" in line:
        Start = int(line.split()[0])
        Stop = int(line.split()[2])
        length = Stop-Start
        Lengths.append(length)
        if length <=10:
            m10.append(length)
        elif length <= 50:
            m50.append(length)
        elif length <= 100:
            m100.append(length)
        elif length <= 200:
            m200.append(length)
        elif length <= 350:
            m350.append(length)
        elif length <= 500:
            m500.append(length)
        elif length <= 650:
            m650.append(length)
        elif length <= 800:
            m800.append(length)
        elif length <= 1000:
            m1000.append(length)



print sum(Lengths) / float(len(Lengths))
print len(m10)
print len(m50)
print len(m100)
print len(m200)
print len(m350)
print len(m500)
print len(m650)
print len(m800)
print len(m1000)












