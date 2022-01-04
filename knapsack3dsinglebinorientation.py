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
    widths = [3,1,2,1,2]
    lengths = [3,1,1,3,2]
    heights=[1,1,2,2,2]
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
                cube=((x >= cube_origins['x'][i]) & (x <= cube_origins['x'][i]+dimension_values_list[0][i]-1) & 
                (y >= cube_origins['y'][i]) & (y <= cube_origins['y'][i]+dimension_values_list[1][i]-1) & 
                (z >= cube_origins['z'][i]) & (z <= cube_origins['z'][i]+dimension_values_list[2][i]-1))
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
    print('creating orientations...')
    orientation_dict=orientations(data['widths'],data['lengths'],data['heights'])
    print('orientations: ', len(orientation_dict))
    max_bin_value=0
    for j in data['bins']:
        max_bin_value+=data['bin_max_width'][j]*data['bin_max_length'][j]*data['bin_max_height'][j]
    max_item_value=0
    for i in data['items']:
        max_item_value+=data['widths'][i]*data['lengths'][i]*data['heights'][i]
    buffer={}
    for h in orientation_dict:
        print('orientation: ',h)
        #print('initializing variables..')
        # Create the mip solver with the SCIP backend.
        solver = pywraplp.Solver.CreateSolver('SCIP')
        
        x={}
        y={}    
        z={}
        ind = {}
        for i in data['items']:
            for j in data['bins']:
                x[(i,j)]=solver.IntVar(0, data['bin_max_width'][j],'x'+str(i)+str(j))
                y[(i,j)]=solver.IntVar(0, data['bin_max_length'][j],'y'+str(i)+str(j))
                z[(i,j)]=solver.IntVar(0, data['bin_max_height'][j],'z'+str(i)+str(j))
                ind[(i, j)] = solver.IntVar(0, 1, 'ind_%i_%i' % (i, j))
        
        x_bool={}
        y_bool={}
        z_bool={}
        M=int(1e6)
        for i in data['items']:
            for j in data['bins']:
                for k in data['items']:
                    if i>=k:
                        continue
                    else:
                        x_bool[(i,j,k)] = solver.IntVar(0,1,"x_bool_"+str(i)+str(j)+str(k))
                        y_bool[(i,j,k)] = solver.IntVar(0,1,"y_bool_"+str(i)+str(j)+str(k))
                        z_bool[(i,j,k)] = solver.IntVar(0,1,"z_bool_"+str(i)+str(j)+str(k))
                        x_bool[(k,j,i)] = solver.IntVar(0,1,"x_bool_"+str(k)+str(j)+str(i))
                        y_bool[(k,j,i)] = solver.IntVar(0,1,"y_bool_"+str(k)+str(j)+str(i))
                        z_bool[(k,j,i)] = solver.IntVar(0,1,"z_bool_"+str(k)+str(j)+str(i))
                        solver.Add(x[(i,j)]+data['widths'][i]<=x[(k,j)]+M*x_bool[(i,j,k)])
                        solver.Add(y[(i,j)]+data['lengths'][i]<=y[(k,j)]+M*y_bool[(i,j,k)])
                        solver.Add(z[(i,j)]+data['heights'][i]<=z[(k,j)]+M*z_bool[(i,j,k)])
                        solver.Add(x[(k,j)]+data['widths'][k]<=x[(i,j)]+M*x_bool[(k,j,i)])
                        solver.Add(y[(k,j)]+data['lengths'][k]<=y[(i,j)]+M*y_bool[(k,j,i)])
                        solver.Add(z[(k,j)]+data['heights'][k]<=z[(i,j)]+M*z_bool[(k,j,i)])
                        tmp=[]
                        tmp.append(x_bool[(i,j,k)])
                        tmp.append(y_bool[(i,j,k)])
                        tmp.append(z_bool[(i,j,k)])
                        tmp.append(x_bool[(k,j,i)])
                        tmp.append(y_bool[(k,j,i)])
                        tmp.append(z_bool[(k,j,i)])
                        solver.Add(sum(tmp)<=5)
            #print('max coordinate constraints...')
            for j in data['bins']:
                solver.Add(x[(i,j)]+orientation_dict[h][0][i]<=data['bin_max_width'][j]+M*(1-ind[(i,j)]))
                solver.Add(y[(i,j)]+orientation_dict[h][1][i]<=data['bin_max_length'][j]+M*(1-ind[(i,j)]))
                solver.Add(z[(i,j)]+orientation_dict[h][2][i]<=data['bin_max_height'][j]+M*(1-ind[(i,j)]))
        
    
        #print('Adding bin indicator constrains...')
        for i in data['items']:
            solver.Add(sum(ind[(i, j)] for j in data['bins']) <= 1)
            
        #print('Adding Volume Constraints...')
        for j in data['bins']:    
            solver.Add(sum(ind[(i,j)]*data['values'][i] for i in data['items'])<=
                       data['bin_max_width'][j]*data['bin_max_length'][j]*data['bin_max_height'][j])
        
    
        # Objective
        objective = solver.Objective()
    
        C=0.01
        objective_terms = []
        for i in data['items']:
            for j in data['bins']:
                objective_terms.append((len(data['bins'])-j)*data['values'][i]*ind[(i,j)])
                objective_terms.append(-C*z[(i,j)]*data['widths'][i]*data['lengths'][i])
        solver.Maximize(solver.Sum(objective_terms))
    
        status = solver.Solve()
        
        if status == pywraplp.Solver.OPTIMAL:
            buffer[h]={}
            buffer[h]['Objective Value']=objective.Value()
            print('Total objective value:', objective.Value())
            total_width = 0
            total_length = 0
            total_height = 0
            total_value = 0
            cube_origins={}
            cube_origins['x']=[]
            cube_origins['y']=[]
            cube_origins['z']=[]
            for j in data['bins']:
                buffer[h]['Bin'+str(j)]={}
                bin_width = 0
                bin_length=0
                bin_height=0
                bin_value = 0
                #print('Bin ', j, '\n')
                for i in data['items']:
                    buffer[h]['Bin'+str(j)]['Item'+str(i)]={}
                    if ind[(i,j)].solution_value() > 0:
                        buffer[h]['Bin'+str(j)]['Item'+str(i)]['Indicator']=1
#                        print('Item', i, '- width:', orientation_dict[h][0][i], ' - length:', orientation_dict[h][1][i],' - height:',
#                              orientation_dict[h][2][i],' value:',data['values'][i])
                        
                        bin_width += orientation_dict[h][0][i]
                        bin_length += orientation_dict[h][1][i]
                        bin_height += orientation_dict[h][2][i]
                        bin_value += data['values'][i]
#                        print(x[(i,j)].solution_value())
#                        print(y[(i,j)].solution_value())
#                        print(z[(i,j)].solution_value())
                        cube_origins['x'].append(x[(i,j)].solution_value())
                        cube_origins['y'].append(y[(i,j)].solution_value())
                        cube_origins['z'].append(z[(i,j)].solution_value())
                    else:
                        buffer[h]['Bin'+str(j)]['Item'+str(i)]['Indicator']=0
                        cube_origins['x'].append(-1)
                        cube_origins['y'].append(-1)
                        cube_origins['z'].append(-1)
#                        print('Item', i, '- width:', orientation_dict[h][0][i], ' - length:', orientation_dict[h][1][i],' - height:',
#                              orientation_dict[h][2][i],' value:',data['values'][i], 'not in bin')
                        
                    buffer[h]['Bin'+str(j)]['Item'+str(i)]['width']=orientation_dict[h][0][i]
                    buffer[h]['Bin'+str(j)]['Item'+str(i)]['length']=orientation_dict[h][1][i]
                    buffer[h]['Bin'+str(j)]['Item'+str(i)]['height']=orientation_dict[h][2][i]
                    buffer[h]['Bin'+str(j)]['Item'+str(i)]['values']=data['values'][i]
#                print('Packed bin value:', bin_value)
                total_width += bin_width
                total_length += bin_length
                total_height += bin_height
                total_value += bin_value
            #make_boxes(data,cube_origins)
            print('Total Packed Value: ', total_value)
        else:
            print(status)
            print('The problem does not have an optimal solution.')
        
        buffer[h]['origin']=cube_origins
        if((total_value==max_bin_value) | (total_value==max_item_value)):
            print('max value achieved: ', total_value)
            print('exiting...')
            break
        
    values=[buffer[b]['Objective Value'] for b in range(0,len(buffer))]
    idx=np.argmax(values)
    print('Solution at index: ', idx)
    print(buffer[idx])
    make_boxes(data,buffer[idx]['origin'],orientation_dict[idx])
# %%
if __name__ == '__main__':
    main()
