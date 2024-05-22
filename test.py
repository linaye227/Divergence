f = open('input.txt','r')
out = open('output.grml', 'w')
lines = f.readlines()
id = 0
out.write('<?xml version="1.0" encoding="UTF-8"?>\n<model formalismUrl="http://formalisms.cosyverif.org/sptgd-net.fml" xmlns="http://cosyverif.org/ns/model">\n')
for line in lines:
    s = line.strip()
    l = [x for x in s.split()]
    if (l[0]=="C"):
        for j in range(1,len(l)):
            out.write('<node id="{}" nodeType="place"> \n\t<attribute name="name">P{}</attribute><attribute name="marking"><attribute name="expr"><attribute name="numValue">0</attribute></attribute></attribute>\n</node>\n'.format(id,j-1))
            id+=1
    elif (l[0]=="S"):
        for j in range(1,len(l)):
            out.write('<node id="{}" nodeType="place"> \n\t<attribute name="name">S{}</attribute><attribute name="marking"><attribute name="expr"><attribute name="numValue">0</attribute></attribute></attribute>\n</node>\n'.format(id,j-1))
            id+=1
f.close()
out.close()
