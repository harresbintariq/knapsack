# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 10:39:47 2021

@author: harres.tariq
"""
from ortools.linear_solver import pywraplp
#import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
# %%
def create_data_model():
    """Create the data for the example."""
    data = {}
    widths = [1,2,1,2]
    lengths = [2,3,1,2]
    values = [widths[i]*lengths[i] for i in range(len(widths))]#[1 for i in range(len(widths))]
    data['widths'] = widths
    data['lengths']=lengths
    data['values'] = values
    data['items'] = list(range(len(widths)))
    data['num_items'] = len(widths)
    num_bins = 1
    data['bins'] = list(range(num_bins))
    data['bin_max_width'] = [3]
    data['bin_max_length'] = [3]
    return data
# %%
def make_rectangles(data,cube_origins,binj):
    clr=['firebrick','yellow','green', 'red','blue','magenta','pink','cyan','khaki','gold','rosybrown','lightcoral','tan','olive']
    fig, ax = plt.subplots()
    ax.set_xlim([0,data['bin_max_width'][binj]])
    ax.set_ylim([0,data['bin_max_length'][binj]])
    for i in data['items']:
        try:
            if(cube_origins[binj]['x'][i]>=0):
                rect= patches.Rectangle((cube_origins[binj]['x'][i], cube_origins[binj]['y'][i]), data['widths'][i], data['lengths'][i], linewidth=5, facecolor=clr[i])
                ax.add_patch(rect)
        except:
            continue
        
    plt.show()
# %%
def main():
    data = create_data_model()

    # Create the mip solver with the SCIP backend.
    solver = pywraplp.Solver.CreateSolver('SCIP')

    # Variables
    # x[i, j] = 1 if item i is packed in bin j.
    x = {}
    for i in data['items']:
        x[i]=solver.IntVar(0, data['bin_max_width'][0],'x'+str(i))
        
    y = {}
    for i in data['items']:
        y[i]=solver.IntVar(0, data['bin_max_length'][0],'y'+str(i))
        
    ind = {}
    for i in data['items']:
        for j in data['bins']:
            ind[(i, j)] = solver.IntVar(0, 1, 'ind_%i_%i' % (i, j))
    
    x_bool={}
    y_bool={}
    M=int(1e6)
    for i in data['items']:
        for j in data['bins']:
            for k in data['items']:
                if i>=k:
                    continue
                else:
                    x_bool[(i,k)] = solver.IntVar(0,1,"x_bool_"+str(i)+str(k))
                    y_bool[(i,k)] = solver.IntVar(0,1,"y_bool_"+str(i)+str(k))  
                    x_bool[(k,i)] = solver.IntVar(0,1,"x_bool_"+str(k)+str(i))
                    y_bool[(k,i)] = solver.IntVar(0,1,"y_bool_"+str(k)+str(i))
                    solver.Add(x[i]+data['widths'][i]<=x[k]+M*x_bool[(i,k)])
                    solver.Add(y[i]+data['lengths'][i]<=y[k]+M*y_bool[(i,k)])
                    solver.Add(x[k]+data['widths'][k]<=x[i]+M*x_bool[(k,i)])
                    solver.Add(y[k]+data['lengths'][k]<=y[i]+M*y_bool[(k,i)])
                    tmp=[]
                    tmp.append(x_bool[(i,k)])
                    tmp.append(y_bool[(i,k)])
                    tmp.append(x_bool[(k,i)])
                    tmp.append(y_bool[(k,i)])
                    solver.Add(sum(tmp)<=3)
            

        for j in data['bins']:
            solver.Add(x[i]+data['widths'][i]<=data['bin_max_width'][0]+M*(1-ind[(i,j)]))
            solver.Add(y[i]+data['lengths'][i]<=data['bin_max_length'][0]+M*(1-ind[(i,j)]))
    

    
    for i in data['items']:
        solver.Add(sum(ind[(i, j)] for j in data['bins']) <= 1)
    for j in data['bins']:    
        solver.Add(sum(ind[(i,j)]*data['values'][i] for i in data['items'])<=
               data['bin_max_width'][0]*data['bin_max_length'][0])
    

    # Objective
    objective = solver.Objective()

    for i in data['items']:
        for j in data['bins']:
            objective.SetCoefficient(ind[(i,j)], data['values'][i])
    objective.SetMaximization()

    status = solver.Solve()

    if status == pywraplp.Solver.OPTIMAL:
        print('Total packed value:', objective.Value())
        total_width = 0
        total_length = 0
        cube_origins={}
        for j in data['bins']:
            bin_width = 0
            bin_length=0
            bin_value = 0
            cube_origins[j]={}
            cube_origins[j]['x']=[]
            cube_origins[j]['y']=[]
            print('Bin ', j, '\n')
            for i in data['items']:
                if ind[(i,j)].solution_value() > 0:
                    print('Item', i, '- width:', data['widths'][i], ' - length:', data['lengths'][i],
                          ' value:',data['values'][i])
                    bin_width += data['widths'][i]
                    bin_length += data['lengths'][i]
                    bin_value += data['values'][i]
                    print(x[i].solution_value())
                    print(y[i].solution_value())
                    cube_origins[j]['x'].append(x[i].solution_value())
                    cube_origins[j]['y'].append(y[i].solution_value())
                    for k in data['items']:
                        if i>=k:
                            continue
                        else:
                            print('i,k: ',i,k)
                            print(x_bool[(i,k)].solution_value())
                            print(y_bool[(i,k)].solution_value())
                            print(x_bool[(k,i)].solution_value())
                            print(y_bool[(k,i)].solution_value())
                else:
                    cube_origins[j]['x'].append(-1)
                    cube_origins[j]['y'].append(-1)
                    print('Item', i, '- width:', data['widths'][i], ' - length:', data['lengths'][i],' value:',data['values'][i],'not in bin')
            print('Packed bin width:', bin_width)
            print('Packed bin length:', bin_length)
            print('Packed bin value:', bin_value)
            print()
            total_width += bin_width
            total_length += bin_length
            print(cube_origins)
            make_rectangles(data,cube_origins,j)
        print('Total packed width:', total_width)
        print('Total packed length:', total_length)
    else:
        print(status)
        print('The problem does not have an optimal solution.')
# %%
if __name__ == '__main__':
    main()

