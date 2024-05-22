f = open('input.txt','r')
out = open('output.grml', 'w')
lines = f.readlines()
id = 0
out.write('<?xml version="1.0" encoding="UTF-8"?>\n<model formalismUrl="http://formalisms.cosyverif.org/sptgd-net.fml" xmlns="http://cosyverif.org/ns/model">\n')
transitioncount = 0
dc={}
ds={}
for line in lines:
    s = line.strip()
    l = [x for x in s.split()]
    if (l[0]=="C"):
        for j in range(1,len(l)):
            out.write('<node id="{}" nodeType="place"> \n\t<attribute name="name">P{}</attribute><attribute name="marking"><attribute name="expr"><attribute name="numValue">0</attribute></attribute></attribute>\n</node>\n'.format(id,j-1))
            dc[j-1] = id
            id+=1
    elif (l[0]=="S"):
        for j in range(1,len(l)):
            out.write('<node id="{}" nodeType="place"> \n\t<attribute name="name">S{}</attribute><attribute name="marking"><attribute name="expr"><attribute name="numValue">0</attribute></attribute></attribute>\n</node>\n'.format(id,j-1))
            ds[j-1] = id
            id+=1
    elif (l[0]=='T'):
        out.write('<node id="{}" nodeType="transition"> \n\t<attribute name="name">T{}</attribute>\n\t<attribute name="priority"><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute>\n\t<attribute name="weight"><attribute name="expr"><attribute name="numValue">{}</attribute></attribute></attribute>\n\t<attribute name="isHorizontal"><attribute name="expr"><attribute name="boolValue">false</attribute></attribute></attribute>\n\t<attribute name="distribution"><attribute name="type">EXPONENTIAL</attribute><attribute name="param"><attribute name="number">0</attribute><attribute name="expr"><attribute name="numValue">1.000000</attribute></attribute></attribute></attribute>\n</node>\n'.format(id,transitioncount,l[5]))
        id+=1
        out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute>\n</arc>\n'.format(id,ds[int(l[1])],id-1))
        id+=1
        out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute>\n</arc>\n'.format(id,dc[int(l[2][1])],id-2))
        id+=1
        out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute>\n</arc>\n'.format(id,id-3,ds[int(l[4])]))
        id+=1
        if (l[3][1]!='-'): out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute>\n</arc>\n'.format(id,id-4,dc[int(l[3][1])]))
        transitioncount+=1
        id+=1
f.close()
out.close()
