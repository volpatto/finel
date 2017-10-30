# -*- coding: utf-8 -*-
import numpy as np
import sys

def fread(file_name, var_str, var_kind):
    '''
    Esta função lê o arquivo file_name (string), busca a string var_str (string) e lê o valor associado e o retorna. A variável var_kind
    especifica o tipo da variável a ser lida.

    Exemplo:
    No arquivo, tem-se:

    var_qualquer{
        23.5
    }

    Fazendo:
    var = fread('nomearquivo.txt','var_qualquer','var_tipo'), var recebe 23.5.

    Autor: Diego T. Volpatto
    E-mail: dtvolpatto@gmail.com
    '''

    f = open(file_name, 'r')
    read_data = f.readlines()
    var = []
    if (var_kind == 'float'):
        for line in range(len(read_data)):
            if ('#' in read_data[line]): 
                continue
            if (var_str in read_data[line]):
                #var = float(read_data[line+1])
                while (read_data[line+1] != '}\n'):
                    #var.append(float(read_data[line+1]))
                    temp_data = read_data[line+1].split()
                    if len(temp_data) > 1:
                        write_data = []
                        for i in range(len(temp_data)):
                            write_data.append(float(temp_data[i]))
                        var.append(write_data)
                    else:
                        var.append(float(read_data[line+1]))
                    line = line + 1
                return var
    if (var_kind == 'int'):
        for line in range(len(read_data)):
            if ('#' in read_data[line]): 
                continue
            if (var_str in read_data[line]):
                #var = int(read_data[line+1])
                while (read_data[line+1] != '}\n'):
                    #var.append(int(read_data[line+1]))
                    temp_data = read_data[line+1].split()
                    if len(temp_data) > 1:
                        write_data = []
                        for i in range(len(temp_data)):
                            write_data.append(int(temp_data[i]))
                        var.append(write_data)
                    else:
                        var.append(int(read_data[line+1]))
                    line = line + 1
                return var
    if (var_kind == 'bool'):
        for line in range(len(read_data)):
            if ('#' in read_data[line]): 
                continue
            if (var_str in read_data[line]):
                #var = bool(read_data[line+1])
                while (read_data[line+1] != '}\n'):
                    var.append(bool(read_data[line+1]))
                    line = line + 1
                return var
    if (var_kind == 'str'):
        for line in range(len(read_data)):
            if ('#' in read_data[line]): 
                continue
            if (var_str in read_data[line]):
                #var = bool(read_data[line+1])
                while (read_data[line+1] != '}\n'):
                    var.append(str(read_data[line+1]).replace('\n',''))
                    line = line + 1
                return var
def getboundary(x,y,bc):
    for i in range(len(bc)):
        if (i!=(len(bc)-1)):
            if (bc[i][1]<=bc[i+1][1] and bc[i][2]<=bc[i+1][2]):
                if ((x>=bc[i][1] and x<=bc[i+1][1]) and (y>=bc[i][2] and y<=bc[i+1][2])):
                    return int(bc[i][0])
            if (bc[i][1]>=bc[i+1][1] and bc[i][2]<=bc[i+1][2]):
                if ((x<=bc[i][1] and x>=bc[i+1][1]) and (y>=bc[i][2] and y<=bc[i+1][2])):
                    return int(bc[i][0])
            if (bc[i][1]<=bc[i+1][1] and bc[i][2]>=bc[i+1][2]):
                if ((x>=bc[i][1] and x<=bc[i+1][1]) and (y<=bc[i][2] and y>=bc[i+1][2])):
                    return int(bc[i][0])
            if (bc[i][1]>=bc[i+1][1] and bc[i][2]>=bc[i+1][2]):
                if ((x<=bc[i][1] and x>=bc[i+1][1]) and (y<=bc[i][2] and y>=bc[i+1][2])):
                    return int(bc[i][0])
        else:
            if (bc[i][1]<=bc[0][1] and bc[i][2]<=bc[0][2]):
                if ((x>=bc[i][1] and x<=bc[0][1]) and (y>=bc[i][2] and y<=bc[0][2])):
                    return int(bc[i][0])
            if (bc[i][1]>bc[0][1] and bc[i][2]<=bc[0][2]):
                if ((x<=bc[i][1] and x>=bc[0][1]) and (y>=bc[i][2] and y<=bc[0][2])):
                    return int(bc[i][0])
            if (bc[i][1]<=bc[0][1] and bc[i][2]>=bc[0][2]):
                if ((x>=bc[i][1] and x<=bc[0][1]) and (y<=bc[i][2] and y>=bc[0][2])):
                    return int(bc[i][0])
            if (bc[i][1]>=bc[0][1] and bc[i][2]>=bc[0][2]):
                if ((x<=bc[i][1] and x>=bc[0][1]) and (y<=bc[i][2] and y>=bc[0][2])):
                    return int(bc[i][0])
    return 0

filename = "usquare.d" 
xdata = fread(filename, "xrange", "float")
ydata = fread(filename, "yrange", "float")
bcdata = fread(filename, "boundaries", "float")
x_pts = np.linspace(xdata[0][0],xdata[0][1],xdata[0][2]+1)
y_pts = np.linspace(ydata[0][0],ydata[0][1],ydata[0][2]+1)
nnodes = len(x_pts)*len(y_pts)
Nx = int(xdata[0][2])
Ny = int(ydata[0][2])
nelem = Nx*Ny
mat = 1
#print bcdata, bcdata[1][1]
#xx, yy = np.meshgrid(x_pts,y_pts)

# Node file
node = open('usquare.n','wt')
print >>node, nnodes
#for i in range(nelem):
for j in range(Ny+1):
    for i in range(Nx+1):
        e = i + j*(Nx+1)
        print >>node, str(e)+':', x_pts[i], y_pts[j], getboundary(x_pts[i],y_pts[j],bcdata)
node.close()

# Element file
elem = open('usquare.e','wt')
print >>elem, nelem
#for i in range(nelem-1):
for j in range(Ny+1):
    for i in range(Nx+1):
        e = i + j*(Nx-1)
        print >>elem, str(e)+':', e, e+1, e+Nx+1, e+Nx,-1, -1, -1, -1, -1, -1, -1, -1, 1
        #print >>elem, str(e)+':', e+1, e+Nx+2, e+Nx+1, e, -1, -1, -1, -1, -1, -1, -1, -1, 1
        #print >>elem, str(e)+':', e, e+1, e+Nx+1, e+Nx+2,-1, -1, -1, -1, -1, -1, -1, -1, 1
elem.close()
sys.exit()
