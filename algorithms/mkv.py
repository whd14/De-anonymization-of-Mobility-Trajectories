import os
from math import *
import argparse
from numpy import *
from collections import Counter


def get_distance(lon1, lat1, lon2, lat2):
    # return distance (kilometer) of two GPS locations
    lon1 = float(lon1)
    lat1 = float(lat1)
    lon2 = float(lon2)
    lat2 = float(lat2)
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    r = 6378.137
    return c * r


def mkv_weight(loc_trace, isp_traces):
    try:
        loc_points = dict(zip([k.strip().split(",")[0] for k in loc_trace.split(";")[1].split("|") if len(k.strip()) > 0],
                              [k.strip().split(",")[1] for k in loc_trace.split(";")[1].split("|") if len(k.strip()) > 0]))
        loc_rg_dict = dict()
        for idx in loc_points.keys():
            min_dis = 1.001
            min_rg = ''
            loc_loc = loc_points[idx]
            for rg in region_c_dict.keys():
                loc = region_c_dict[rg]
                dis = get_distance(float(loc.split("_")[0]), float(loc.split("_")[1]), float(loc_loc.split("_")[1]),
                                   float(loc_loc.split("_")[0]))
                if dis < min_dis:
                    min_dis = dis
                    min_rg = rg
            if region_c_dict.has_key(min_rg):
                loc_rg_dict[idx]=region_c_dict[min_rg]
    except IndexError:
        print 'loc_trace', loc_trace
    weight_list = []
    for isp_trace in isp_traces:
        isp_id = isp_trace.split(";")[0]
        weight=mkv_single(loc_rg_dict,isp_trace)
        weight_list.append((isp_id,weight))
    return weight_list

def mkv_single(loc_rg_dict, isp_trace):
    p_zero=10**(-100) #Smooth parameter to eliminate zero probability
    try:
        isp_points = dict(zip([k.strip().split(",")[0] for k in isp_trace.split(";")[1].split("|") if len(k.strip()) > 0 and len(k.split(",")) >= 2],
                              [region_c_dict[k.strip().split(",")[1]] for k in isp_trace.split(";")[1].split("|") if len(k.strip()) > 0 and len(k.split(",")) >= 2]))
    except IndexError:
        print 'loc_trace', loc_trace
        return float(-Inf)
        
    RegA=list(set(isp_points.values()+loc_rg_dict.values()))
    Nrg=len(RegA)
    RegMap=dict(zip(RegA,range(0,Nrg)))
    
    Ezero=float(0.5) #Parameters to eliminate zero transition and marginal probability
    Tzero=Ezero/Nrg

    EdgeISP=zeros(Nrg+1) #Marginal distribution of ISP trajectory
    NlogISP=len(isp_points.values())
    LogHistISP=Counter(isp_points.values())

    for rg in RegA:
        if LogHistISP.has_key(rg):
            EdgeISP[RegMap[rg]]=float(Ezero+LogHistISP[rg])/(NlogISP+Ezero*(Nrg+1))
        else:
            EdgeISP[RegMap[rg]]=float(Ezero)/(NlogISP+Ezero*(Nrg+1))
    EdgeISP[Nrg]=float(Ezero)/(NlogISP+Ezero*(Nrg+1))
    TransISP = zeros([Nrg+1]*2)
    TimeSlotsISP = [int(x) for x in isp_points.keys()]
    TimeSlotsISP.sort()
    for n in range(0,len(TimeSlotsISP)-1):
        rstart=RegMap[isp_points[str(TimeSlotsISP[n])]]
        dstart=RegMap[isp_points[str(TimeSlotsISP[n+1])]]
        TransISP[rstart,dstart]=TransISP[rstart,dstart]+1
    for sint in range(0,Nrg):
        TransISP[sint,:]=(TransISP[sint,:]+Tzero)/(sum(TransISP[sint,:])+Tzero*(Nrg+1))
    TransISP[Nrg,:]=EdgeISP #Transition probability matrix of ISP trajectory
    
    EdgeLOC=zeros(Nrg+1) #Marginal distribution of external trajectory
    NlogLOC=len(loc_rg_dict.values())
    LogHistLOC=Counter(loc_rg_dict.values())
    for rg in RegA:
        if LogHistLOC.has_key(rg):
            EdgeLOC[RegMap[rg]]=float(Ezero+LogHistLOC[rg])/(NlogLOC+Ezero*(Nrg+1))
        else:
            EdgeLOC[RegMap[rg]]=float(Ezero)/(NlogLOC+Ezero*(Nrg+1))
    EdgeLOC[Nrg]=float(Ezero)/(NlogLOC+Ezero*(Nrg+1))
    TransLOC = zeros([Nrg+1]*2) #Transition probability matrix of external trajectory
    TimeSlotsLOC = [int(x) for x in loc_rg_dict.keys()]
    TimeSlotsLOC.sort()
    for n in range(0,len(TimeSlotsLOC)-1):
        rstart=RegMap[loc_rg_dict[str(TimeSlotsLOC[n])]]
        dstart=RegMap[loc_rg_dict[str(TimeSlotsLOC[n+1])]]
        TransLOC[rstart,dstart]=TransLOC[rstart,dstart]+1
    for sint in range(0,Nrg):
        TransLOC[sint,:]=(TransLOC[sint,:]+Tzero)/(sum(TransLOC[sint,:])+Tzero*(Nrg+1))
    TransLOC[Nrg,:]=EdgeLOC
    weight=multiply(TransISP.T*EdgeISP,TransLOC.T*EdgeLOC).sum()
    return weight


def hit_precision(rank, k):
    # return hit-precision with parameter K=k
    if rank > k:
        rank = k
    return 1 - rank * 1.0 / k


def single_task(locfile_name, input_path, output_path, save_weights=False):
    # deanonymization task of single user
    fr = open(input_path + locfile_name, 'r')
    all_traces = fr.readlines()
    loc_trace = all_traces[0]
    isp_traces = [all_traces[k] for k in range(len(all_traces)) if k >= 1]
    ispUID_score = mkv_weight(loc_trace, isp_traces)
    ispUID_score_sorted = sorted(ispUID_score, key=lambda x: x[1], reverse=True)
    ispUID_sorted = [x[0] for x in ispUID_score_sorted]
    locUID = loc_trace.split(";")[0]
    if locUID in ispUID_sorted:
        rank_true = ispUID_sorted.index(locUID)
    else:
        rank_true = len(ispUID_sorted)
    acc = hit_precision(rank_true, K_hp)
    if save_weights:
        fw = open(output_path + locfile_name, 'w')
        for i in range(len(ispUID_score_sorted)):
            fw.write(str(ispUID_score_sorted[i][0]) + '\t' + str(ispUID_score_sorted[i][1]) + '\n')
        fw.close()
    return locfile_name, acc


def main(input_path='./input/', output_path='./output/', save_weights=False):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    file_names = os.listdir(input_path)
    acc_list = []
    for file_name in file_names:
        acc_list.append(single_task(file_name, input_path, output_path, save_weights))
    with open(output_path + '/result', 'w') as fw:
        for acc in acc_list:
            fw.write(str(acc[0]) + '\t' + str(acc[1]) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', nargs='?', default='./input/', help='Path of input data')
    parser.add_argument('--output_path', nargs='?', default='./output/mkv/', help='Path of output data')
    parser.add_argument('--region_centre', nargs='?', default='./RegionCenters', help='Path of region centres\' location ')
    parser.add_argument('--save_weights', nargs='?', default=False, help='Save weight or not')
    parser.add_argument('--scale', nargs='?', default=0.01, help='Spatial scale', type=float)
    parser.add_argument('--K', nargs='?', default=10, help='K in hit-precision', type=int)

    args = parser.parse_args()
    K_hp = args.K
    region_centre = args.region_centre
    region_c_list = open(region_centre, 'r').readlines()
    region_c_dict = dict(
        zip([k.split("_")[0] for k in region_c_list], [k.split("_")[1] + "_" + k.split("_")[2] for k in region_c_list]))

    centre_map_dict = {}
    for centre in region_c_dict:
        loc = region_c_dict[centre].split('_')
        new_loc = str(floor(float(loc[0]) / args.scale) * args.scale) + '-' + str(floor(float(loc[1]) / args.scale) * args.scale)
        centre_map_dict[centre] = new_loc

    main(input_path=args.input_path, output_path=args.output_path, save_weights=True)
