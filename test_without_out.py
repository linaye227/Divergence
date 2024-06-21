def printGRML_OneTransition(id,name,distribution,para,priority,weight,extra):
    s=""
    s+='  <node id=\"'+ id + '\"  nodeType=\"transition\">\n'
    s+='    <attribute name=\"name\">'+ name +'</attribute>\n'
    s+='    <attribute name=\"distribution\"> <attribute name=\"type\">\n'
    s+=distribution
    s+= '</attribute><attribute name="param"><attribute name="number">0</attribute><attribute name="expr">'
    s+=para
    s+='    </attribute></attribute>\n'
    s+='    </attribute>\n'
    s+='    <attribute name=\"priority\"><attribute name=\"expr\">\n'
    s+='    <attribute name=\"numValue\">'+ priority + '</attribute>\n'
    s+='    </attribute></attribute>\n'
    s+='    <attribute name=\"weight\"><attribute name=\"expr\">\n'
    s+='    '+ '<attribute name="numValue">' + weight + '</attribute>' +  '\n'
    s+='    </attribute></attribute>\n'
    s+='    '+extra+'\n'
    s+='  </node>\n'
    return(s)
def poly(n,l):
    k = 1
    ans = 0
    for i in range(len(l)):
        ans+=k*l[i]
        k*=n
    return ans
minimum = 1e6
f = open('input.txt','r')
out = open('output.grml', 'w')
lines = f.readlines()
id = 0
dw = {}
out.write('<?xml version="1.0" encoding="UTF-8"?>\n<model formalismUrl="http://formalisms.cosyverif.org/sptgd-net.fml" xmlns="http://cosyverif.org/ns/model">\n')
transitioncount = 0
dc={}
ds={}
do = {}
di = {}
L = []
channelcount = 0
out.write('<attribute name="declaration"><attribute name="constants"><attribute name="intConsts"><attribute name="intConst"><attribute name="name">NbClient</attribute><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute><attribute name="intConst"><attribute name="name">NbPosition</attribute><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute></attribute></attribute><attribute name="classDeclaration"><attribute name="name">client</attribute><attribute name="classType"><attribute name="classIntInterval"><attribute name="lowerBound">1</attribute><attribute name="higherBound"><attribute name="name">NbClient</attribute></attribute></attribute></attribute><attribute name="circular">false</attribute></attribute>')
counter = 0
sumw = 1
for line in lines:
    counter+=1
    s = line.strip()
    l = [x for x in s.split()]
    if (l[0]=="C"):
        channelcount = len(l)-1
        for j in range(1,len(l)):
            if (j-1!=0) : 
                out.write('<node id="{}" nodeType="place"> \n\t<attribute name="name">C{}</attribute>\n\t<attribute name="domain"><attribute name="type">CP{}</attribute></attribute>\n\t<attribute name="marking"><attribute name="expr"><attribute name="numValue">0</attribute></attribute></attribute>\n</node>\n'.format(id,j-1,j-1))
                dc[j-1] = id
                id+=1
                out.write('<node id="{}" nodeType="place"> \n\t<attribute name="name">In{}</attribute>\n\t<attribute name="domain"><attribute name="type">position{}</attribute></attribute>\n\t<attribute name="marking"><attribute name="token"><attribute name="occurs"><attribute name="intValue">1</attribute></attribute><attribute name="tokenProfile"><attribute name="expr"><attribute name="enumConst"><attribute name="type">position{}</attribute><attribute name="enumValue">n0</attribute></attribute></attribute></attribute></attribute></attribute>\n</node>\n'.format(id,j-1,j-1,j-1))
                di[j-1] = id
                id+=1
            else :
                for i in range(1,len(l)-1):
                    out.write('<attribute name="classDeclaration"><attribute name="name">position{}</attribute><attribute name="classType"><attribute name="classIntInterval"><attribute name="lowerBound">1</attribute><attribute name="higherBound"><attribute name="name">NbPosition</attribute></attribute></attribute></attribute><attribute name="circular">true</attribute></attribute>'.format(i))
                    out.write('<attribute name="domainDeclaration"><attribute name="name">CP{}</attribute><attribute name="domainType"><attribute name="cartesianProduct"><attribute name="type">client</attribute><attribute name="type">position{}</attribute></attribute></attribute></attribute>'.format(i,i))
                for i in range(1,len(l)-1):
                    out.write('<attribute name="variableDeclaration"><attribute name="name">p{}</attribute><attribute name="type">position{}</attribute></attribute>'.format(i,i))
                out.write('<attribute name="variableDeclaration"><attribute name="name">x</attribute><attribute name="type">client</attribute></attribute></attribute>'.format(i))
                out.write('\n<node id="{}" nodeType="place"> \n\t<attribute name="name">Count</attribute><attribute name="marking"><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute>\n</node>\n'.format(id))
                id+=1
                out.write('<node id="{}" nodeType="place">\n\t<attribute name="name">C{}</attribute>\n\t<attribute name="marking"><attribute name="expr"><attribute name="numValue">0</attribute></attribute></attribute>\n</node>\n'.format(id,j-1))
                dc[j-1] = id
                id+=1
    elif (l[0]=="S"):
        for j in range(1,len(l)):
            dw[j-1] = 0
            if (j==1): out.write('<node id="{}" nodeType="place"> \n\t<attribute name="name">S{}</attribute><attribute name="marking"><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute>\n</node>\n'.format(id,j-1))
            else : out.write('<node id="{}" nodeType="place"> \n\t<attribute name="name">S{}</attribute><attribute name="marking"><attribute name="expr"><attribute name="numValue">0</attribute></attribute></attribute>\n</node>\n'.format(id,j-1))
            ds[j-1] = id
            id+=1
    elif (l[0]=='T'):
        sumw +=int(l[5])
        dw[int(l[1])]+=int(l[5])
        transid = id
        out.write('<node id="{}" nodeType="transition"> \n\t<attribute name="name">T{}</attribute>\n\t<attribute name="priority"><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute>\n\t<attribute name="weight"><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute>\n\t<attribute name="isHorizontal"><attribute name="expr"><attribute name="boolValue">false</attribute></attribute></attribute>\n\t<attribute name="distribution"><attribute name="type">EXPONENTIAL</attribute><attribute name="param"><attribute name="number">0</attribute><attribute name="expr"><attribute name="numValue">{}</attribute></attribute></attribute></attribute>\n</node>\n'.format(id,transitioncount,l[5]))
        id+=1
        out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute>\n</arc>\n'.format(id,ds[int(l[1])],transid))
        id+=1
        out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute>\n</arc>\n'.format(id,transid,ds[int(l[4])]))
        id+=1
        if(l[3][1]=='-'):
            out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute>\n</arc>\n'.format(id,0,transid))
            id+=1
            out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="expr"><attribute name="numValue">0</attribute></attribute></attribute>\n</arc>\n'.format(id,transid,0))
            id+=1
        transitioncount+=1
        tempo = False
        if (l[2][1]!='0'):
            out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="token"><attribute name="occurs"><attribute name="intValue">1</attribute></attribute><attribute name="tokenProfile"><attribute name="expr"><attribute name="name">p{}</attribute></attribute></attribute></attribute></attribute>\n</arc>\n'.format(id,di[int(l[2][1])],transid,int(l[2][1])))
            id+=1
            if ( l[3][1]=='-' or di[int(l[2][1])]!=di[int(l[3][1])]):
                out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="token"><attribute name="occurs"><attribute name="intValue">1</attribute></attribute><attribute name="tokenProfile"><attribute name="expr"><attribute name="function"><attribute name="--"><attribute name="name">p{}</attribute><attribute name="intValue">1</attribute></attribute></attribute></attribute></attribute></attribute></attribute>\n</arc>\n'.format(id,transid,di[int(l[2][1])],int(l[2][1])))
                id+=1
                tempo = True
            else:
                out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="token"><attribute name="occurs"><attribute name="intValue">1</attribute></attribute><attribute name="tokenProfile"><attribute name="expr"><attribute name="name">p{}</attribute></attribute></attribute></attribute></attribute>\n</arc>\n'.format(id,transid,di[int(l[2][1])],int(l[2][1])))
                id+=1
            out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="token"><attribute name="occurs"><attribute name="intValue">1</attribute></attribute><attribute name="tokenProfile"><attribute name="expr"><attribute name="name">x</attribute></attribute><attribute name="expr"><attribute name="name">p{}</attribute></attribute></attribute></attribute></attribute>\n</arc>\n'.format(id,dc[int(l[2][1])],transid,int(l[2][1])))
            id+=1
        else:
            out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute>\n</arc>\n'.format(id,dc[int(l[2][1])],transid))
            id+=1
            tempo = True
        if (l[3][1]!='-'):
            if (tempo):
                out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="token"><attribute name="occurs"><attribute name="intValue">1</attribute></attribute><attribute name="tokenProfile"><attribute name="expr"><attribute name="name">p{}</attribute></attribute></attribute></attribute></attribute>\n</arc>\n'.format(id,di[int(l[3][1])],transid,int(l[3][1])))
                id+=1
                out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="token"><attribute name="occurs"><attribute name="intValue">1</attribute></attribute><attribute name="tokenProfile"><attribute name="expr"><attribute name="function"><attribute name="++"><attribute name="name">p{}</attribute><attribute name="intValue">1</attribute></attribute></attribute></attribute></attribute></attribute></attribute>\n</arc>\n'.format(id,transid,di[int(l[3][1])],int(l[3][1])))
                id+=1
            out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="token"><attribute name="occurs"><attribute name="intValue">1</attribute></attribute><attribute name="tokenProfile"><attribute name="expr"><attribute name="name">x</attribute></attribute><attribute name="expr"><attribute name="name">p{}</attribute></attribute></attribute></attribute></attribute>\n</arc>\n'.format(id,transid,dc[int(l[3][1])],int(l[3][1])))
            id+=1
    else:
        ll = []
        for j in range(1,len(l)):
            ll.append(int(l[j]))
        L.append(ll)
        start = ""
        end = ""
        for k in range(len(l)-2):
            start+= '\t<attribute name="function"><attribute name="+"> \n'
        for j in range (1,len(l)):
            start2 = ""
            end2 = ""
            for k in range(j-1):
                start2+= '\t<attribute name="function"><attribute name="*"> \n'
                end2 = '\t<attribute name="function"><attribute name="-"><attribute name="name">Count</attribute><attribute name="numValue">1</attribute>\n\t</attribute></attribute></attribute></attribute> \n' + end2
            start2+= '\t<attribute name="numValue">' + l[j] + '</attribute>\n'
            start+=start2+end2
            if (j>1):
                start+= '\t</attribute></attribute> \n'
        s = printGRML_OneTransition(str(id),"T"+str(transitioncount),"EXPONENTIAL",start,"0",str(1),"")
        transid = id
        transitioncount+=1
        id+=1
        out.write(s);
        out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute>\n</arc>\n'.format(id,0,transid))
        id+=1
        out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="token"><attribute name="occurs"><attribute name="intValue">1</attribute></attribute><attribute name="tokenProfile"><attribute name="expr"><attribute name="name">x</attribute></attribute></attribute></attribute></attribute>\n</arc>\n'.format(id,transid,1))
        id+=1
        out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="expr"><attribute name="numValue">2</attribute></attribute></attribute>\n</arc>\n'.format(id,transid,0))
        id+=1
        out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute>\n</arc>\n'.format(id,ds[counter-3],transid))
        id+=1
        out.write('<arc id="{}" arcType="arc" source="{}" target="{}">\n\t<attribute name="valuation"><attribute name="expr"><attribute name="numValue">1</attribute></attribute></attribute>\n</arc>\n'.format(id,transid,ds[counter-3]))
        id+=1
for i in range(len(L)):
    l = 0
    r = minimum
    while(r-l>1):
        mid = (l+r)//2
        if (poly(mid,L[i])<sumw): l = mid
        else: r = mid
    minimum = min(minimum,r)
print("N0 =",minimum)
maxiw = 0
for i in dw:
    maxiw = max(maxiw,dw[i])
out.write('\n<node id="{}" nodeType="place"> \n\t<attribute name="name">Wmax</attribute><attribute name="marking"><attribute name="expr"><attribute name="numValue">{}</attribute></attribute></attribute>\n</node>\n'.format(id,maxiw))
id+=1
out.write('\n<node id="{}" nodeType="place"> \n\t<attribute name="name">N0</attribute><attribute name="marking"><attribute name="expr"><attribute name="numValue">{}</attribute></attribute></attribute>\n</node>\n'.format(id,minimum))
id+=1
out.write("</model>\n")
f.close()
out.close()
