# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 10:39:47 2021

@author: harres.tariq
"""
from ortools.linear_solver import pywraplp
import numpy as np
import matplotlib.pyplot as plt
import itertools
# %%
def orientations(widths,lengths,heights):
    tmp=[widths, lengths,heights]
    #list(itertools.permutations([1, 2, 3]))
    a={}
    for i in range(0,len(tmp[0])):
        ls=[]
        ls.append(tmp[0][i])
        ls.append(tmp[1][i])
        ls.append(tmp[2][i])
        a[i]=list(itertools.permutations(ls))
        
    w=np.zeros((6**len(a),len(a)), np.int)
    l=np.zeros((6**len(a),len(a)), np.int)
    h=np.zeros((6**len(a),len(a)), np.int)
    for i in range(0,len(a)):
        threshold=int(6**len(a)/6**(i+1))
        k=0
        for j in range(0,6**(i+1)):
            #print(i,j,k)
            w[j*threshold:(j+1)*threshold,i]=a[i][k][0]
            l[j*threshold:(j+1)*threshold,i]=a[i][k][1]
            h[j*threshold:(j+1)*threshold,i]=a[i][k][2]
            if k==5:
                k=0
            else:
                k+=1
            
    d={}
    for i in range(0, w.shape[0]):      
        d[i]=[list(w[i,:]), list(l[i,:]), list(h[i,:])]
                
            
        
    return d
# %%
def create_data_model():
    """Create the data for the example."""
    data = {}
    widths = [3,2,1,1,2]
    lengths = [3,2,1,2,2]
    heights=[1,1,1,1,2]
    values = [widths[i]*lengths[i]*heights[i] for i in range(len(widths))]#[1 for i in range(len(widths))]
    data['widths'] = widths
    data['lengths']=lengths
    data['heights']=heights
    data['values'] = values
    data['items'] = list(range(len(widths)))
    data['num_items'] = len(widths)
    num_bins = 1
    data['bins'] = list(range(num_bins))
    data['bin_max_width'] = [3]
    data['bin_max_length'] = [3]
    data['bin_max_height']=[3]
    return data
# %%
def make_boxes(data,cube_origins,dimension_values_list):
    x,y,z=np.indices((data['bin_max_width'][0],data['bin_max_length'][0],data['bin_max_height'][0]))
    voxels=(np.zeros((data['bin_max_width'][0],data['bin_max_length'][0],data['bin_max_height'][0]))==True)
    colors = np.empty(voxels.shape, dtype=object)
    #n=len(cube_origins)
    #clr = cm.rainbow(np.linspace(0, 1, n))
    clr=['firebrick','yellow','green', 'red','blue','magenta','pink','cyan','khaki','gold','rosybrown','lightcoral','tan','olive']
       
    for i in data['items']:
        try:
            if(cube_origins['x'][i]>=0):
                cube=(x >= cube_origins['x'][i]) & (x <= cube_origins['x'][i]+dimension_values_list[0][i]-1) & (y >= cube_origins['y'][i]) & (y <= cube_origins['y'][i]+dimension_values_list[1][i]-1) & (z >= cube_origins['z'][i]) & (z <= cube_origins['z'][i]+dimension_values_list[2][i]-1)
                colors[cube]=clr[i]
                voxels |= cube
        except:
            continue

    ax = plt.figure().add_subplot(projection='3d')
    ax.voxels(voxels, facecolors=colors, edgecolor='k')

    plt.show()
# %%
def main():
    data = create_data_model()
    
    # Create the mip solver with the SCIP backend.
    solver = pywraplp.Solver.CreateSolver('SCIP')
    print('creating orientations...')
    orientation_dict=orientations(data['widths'],data['lengths'],data['heights'])
    print('orientations: ', len(orientation_dict))
    print('initializing variables..')
    ori_ind={}
    x={}
    y={}    
    z={}
    ind = {}
    for h in orientation_dict:
        print(h)
        ori_ind[h]=solver.IntVar(0,1,'ori_ind'+str(h))
        for i in data['items']:
            for j in data['bins']:
                x[(h,i,j)]=solver.IntVar(0, data['bin_max_width'][j],'x'+str(h)+str(i)+str(j))
                y[(h,i,j)]=solver.IntVar(0, data['bin_max_length'][j],'y'+str(h)+str(i)+str(j))
                z[(h,i,j)]=solver.IntVar(0, data['bin_max_height'][j],'z'+str(h)+str(i)+str(j))
                ind[(h,i,j)] = solver.IntVar(0, 1, 'ind_%i_%i_%i' % (h,i,j))
            
        x_bool={}
        y_bool={}
        z_bool={}
        M=int(1e6)
        print('creating coordinate constraints...')
        for i in range(0,len(data['items'])):
            for j in data['bins']:
                for k in range(i+1, len(data['items'])):
                    print(i,j,k)
                    if i==k:
                        continue
                    else:
                        x_bool[(h,i,j)] = solver.IntVar(0,1,"x_bool_"+str(h)+str(i)+str(j))   
                        y_bool[(h,i,j)] = solver.IntVar(0,1,"y_bool_"+str(h)+str(i)+str(j))
                        z_bool[(h,i,j)] = solver.IntVar(0,1,"z_bool_"+str(h)+str(i)+str(j))
                        solver.Add(x[(h,i,j)]+orientation_dict[h][0][i]<=x[(h,k,j)]+M*(3-x_bool[(h,i,j)]-ind[(h,i,j)]-ori_ind[h]))
                        solver.Add(y[(h,i,j)]+orientation_dict[h][1][i]<=y[(h,k,j)]+M*(3-y_bool[(h,i,j)]-ind[(h,i,j)]-ori_ind[h]))
                        solver.Add(z[(h,i,j)]+orientation_dict[h][2][i]<=z[(h,k,j)]+M*(3-z_bool[(h,i,j)]-ind[(h,i,j)]-ori_ind[h]))
                        tmp=[] 
                        tmp.append(x_bool[(h,i,j)])
                        tmp.append(y_bool[(h,i,j)])
                        tmp.append(z_bool[(h,i,j)])
                        solver.Add(sum(tmp)-1*ori_ind[h]==0)#
            print('creating coordinate max constraints...')
            for j in data['bins']:
                solver.Add(x[(h,i,j)]+orientation_dict[h][0][i]<=data['bin_max_width'][j]+M*(2-ind[(h,i,j)]-ori_ind[h]))
                solver.Add(y[(h,i,j)]+orientation_dict[h][1][i]<=data['bin_max_length'][j]+M*(2-ind[(h,i,j)]-ori_ind[h]))
                solver.Add(z[(h,i,j)]+orientation_dict[h][2][i]<=data['bin_max_height'][j]+M*(2-ind[(h,i,j)]-ori_ind[h]))
        

        print('bin indicator constraints...')
        for i in data['items']:
            solver.Add(sum(ind[(h,i,j)] for j in data['bins']) -1*ori_ind[h]<= 0)#
        
        print('volume constraint...')
        for j in data['bins']:    
            solver.Add(sum(ind[(h,i,j)]*data['values'][i] for i in data['items'])<=
                       data['bin_max_width'][j]*data['bin_max_length'][j]*data['bin_max_height'][j])
    
    print('orientation constraint...')
    solver.Add(sum(ori_ind[h] for h in orientation_dict)==1)
    

    # Objective
    objective = solver.Objective()
    for h in orientation_dict:
        for i in data['items']:
            for j in data['bins']:
                objective.SetCoefficient(ind[(h,i,j)], data['values'][i])
    objective.SetMaximization()
    print('starting optimization...')
    status = solver.Solve()
    print('results...')
    if status == pywraplp.Solver.OPTIMAL:
        print('Total packed value:', objective.Value())
        total_width = 0
        total_length = 0
        total_height = 0
        cube_origins={}
        cube_origins['x']=[]
        cube_origins['y']=[]
        cube_origins['z']=[]
        for h in orientation_dict:
            if(ori_ind[h].solution_value()>0):
                print('orientation: ',h)
                print(orientation_dict[h])
                for j in data['bins']:
                    bin_width = 0
                    bin_length=0
                    bin_height=0
                    bin_value = 0
                    print('Bin ', j, '\n')
                    for i in data['items']:
                        if ind[(h,i,j)].solution_value() > 0:
                            print('Item', i, '- width:', orientation_dict[h][0][i], ' - length:', orientation_dict[h][1][i],' - height:',
                                  orientation_dict[h][2][i],' value:',data['values'][i])
                            bin_width += orientation_dict[h][0][i]
                            bin_length += orientation_dict[h][1][i]
                            bin_height += orientation_dict[h][2][i]
                            bin_value += data['values'][i]
                            print(x[(h,i,j)].solution_value())
                            print(y[(h,i,j)].solution_value())
                            print(z[(h,i,j)].solution_value())
                            cube_origins['x'].append(x[(h,i,j)].solution_value())
                            cube_origins['y'].append(y[(h,i,j)].solution_value())
                            cube_origins['z'].append(z[(h,i,j)].solution_value())
                        else:
                            cube_origins['x'].append(-1)
                            cube_origins['y'].append(-1)
                            cube_origins['z'].append(-1)
                            print('Item', i, '- width:', orientation_dict[h][0][i], ' - length:', orientation_dict[h][1][i],' - height:',
                                  orientation_dict[h][2][i],' value:',data['values'][i],'not in bin')
        
                    print('Packed bin value:', bin_value)
                    print()
                    total_width += bin_width
                    total_length += bin_length
                    total_height += bin_height
                make_boxes(data,cube_origins,orientation_dict[h])
                break

    else:
        print(status)
        print('The problem does not have an optimal solution.')
# %%
if __name__ == '__main__':
    main()
