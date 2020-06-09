#from __future__ import annotations
from itertools import chain
from typing import Union, List, Tuple, Dict, Set, FrozenSet
from time import time, ctime, strftime
from datetime import timedelta
from csv import writer
from functools import partial
import re

class sortedset(list):
  def append(self, other, other_value):
    done = False
    for i, (obj, value) in enumerate(self):
      if not done and other_value < value:
        self.insert(i, (other, other_value))
        done = True
      if obj is other:
        if done:
            self.pop(i)
        return
    if not done:
      super().append((other, other_value))
          

class Vertex:
    def __init__(self, loc: Tuple[float], color: Tuple[float]):
        self.x = loc[0]
        self.y = loc[1]
        self.z = loc[2]
        self.number = 0
        self.area = 0
        self.classification = Classification_From_Color(color)
        self.touching = set()
        self.patch = None
        self.center = False
        self.highest = False

    def reintialize(self, identifier: int):
        self.ids.append(identfier)
        return self

    def distance_to(self, other) -> float:
        return ((self.x-other.x)**2 + (self.y-other.y)**2 + (self.z-other.z)**2)**.5

    def __str__(self) -> str:
        return str(self.x) + str(self.y) + str(self.z) + str(self.classification)


class Patch:
    def __init__(self, vertex: Vertex, rim: bool):
        self.classification = vertex.classification
        self.distances = []
        self.vertices = {vertex,}
        self.rim = set()
        self.core = set()
        (self.rim if rim else self.core).add(vertex)
        self.touching = set()
        self.surface = None
        self.center = None
        self.highest = False

    def finish(self):
        i = 0
        for vertex in self.rim:
            vertex.number = i
            vertex.patch = self
            i += 1
        for vertex in self.core:
            vertex.number = i
            vertex.patch = self
            i += 1
        self.distances = [[0] * len(self.rim) for _ in range(len(self.vertices))] 

    def is_touching_patch(self, other) -> bool:
        return any(any(touching in other.rim for touching in vertex.touching) for vertex in self.rim)

    def distance_to(self, other) -> float:
        return distances_dict[frozenset((self.center, other.center))]

    def add_patch(self, patch):
        if patch.highest: self.highest = True
        self.vertices |= patch.vertices
        self.rim |= patch.rim
        self.core |= patch.core
    
    def get_center(self) -> Vertex:
        if not self.center and self.rim:
                best = (None, float("inf"))
                for vertex in self.vertices:
                    try:
                        distances = self.distances[vertex.number]
                        mean = sum(distances)/len(distances)
                        absolute_deviation = sum((mean - distance)**2 for distance in distances)
                        mean_absolute_deviation = absolute_deviation**.5/len(distances)
                        if mean_absolute_deviation < best[1]:
                            best = (vertex, mean_absolute_deviation)
                    except Exception as e:
                        print("Distances:", distances, self.distances)
                        print("Len:", len(distances))
                        print("Vertex:", vertex)
                        raise e
                self.center = best[0]
                self.center.center = True
            


    def area(self) -> float:
        return sum(vertex.area for vertex in self.vertices)/3

RBG_dict = {
(0.75, 0.0, 0.75): 'Charged -', 
(0.0, 1.0, 0.0): 'Charged +', 
(1.0, 0.0, 0.0): 'Hydrophilic', 
(0.0, 0.0, 1.0): 'Hydrophobic'
}


def Make_Face(distances_dict: dict, vertex1: Vertex, vertex2: Vertex, vertex3: Vertex):
    vertex1.touching.add(vertex2)
    vertex2.touching.add(vertex1)
    a = vertex1.distance_to(vertex2)
    distances_dict[frozenset((vertex1, vertex2))] = a
    vertex2.touching.add(vertex3)
    vertex3.touching.add(vertex2)
    b = vertex2.distance_to(vertex3)
    distances_dict[frozenset((vertex2, vertex3))] = b 
    vertex3.touching.add(vertex1)
    vertex1.touching.add(vertex3)
    c = vertex3.distance_to(vertex1)
    distances_dict[frozenset((vertex3, vertex1))] = c
    p = (a+b+c)/2
    area = (p*(p-a)*(p-b)*(p-c))**.5
    vertex1.area += area 
    vertex2.area += area 
    vertex3.area += area 


def Classification_From_Color(color: Tuple[float]) -> str:
    if color in RBG_dict:
        return RBG_dict[color]
    best = [3,]
    for key in RBG_dict:
        red = abs(key[0] - color[0])
        green = abs(key[1] - color[1])
        blue = abs(key[2] - color[2])
        score = red + green + blue
        if score <= best[0]:
            best = [score, key]
    return RBG_dict[best[1]]

def Read_Wrl(distances_dict: dict, file: str) -> Set[Vertex]:
    try:
        with open(file, 'r') as wrl:
            mode = 0
            vertices_data = []
            faces_data = []
            color_data = []
            for line in wrl:
                if ']' in line:
                    mode = 0
                if not mode:
                    if 'point [' in line:
                        mode = 1
                    elif 'coordIndex [' in line:
                        mode = 2
                    elif 'color [' in line:
                        mode = 3
                    continue
                if mode == 1:
                    a = re.findall("[-0-9]{1,3}.[0-9]{6}", line)
                    vertices_data.append(tuple(map(float, a)))
                elif mode == 2:
                    a = re.findall("[0-9]{1,}", line)
                    faces_data.append(tuple(map(int, a)))
                elif mode == 3:
                    a = re.findall("[0,1].[0-9]{4}", line)
                    color_data.append(tuple(map(float, a)))
            vertex_dict = {}
            vertices = []
            highest = [None, float("-inf")]
            for loc, color in zip(vertices_data, color_data):
                if loc in vertex_dict:
                    vertices.append(vertices[vertex_dict[loc]])
                    continue
                vertex_dict[loc] = len(vertices)
                vertex = Vertex(loc, color)
                vertices.append(vertex)
                if loc[1] > highest[1]:
                    highest = [vertex, loc[1]]
            highest[0].highest = True
            for face in faces_data:
                Make_Face(distances_dict, vertices[face[0]], vertices[face[1]], vertices[face[2]]) 
            vertices = set(vertices)
            vertices.discard(None)
            return vertices
    except Exception:
        raise 


def Get_Patches(vertices: Set[Vertex]) -> Tuple[Set[Patch], int]:
    patches = set()
    for vertex in vertices:
        patchable = set()
        rim = False
        for other in vertex.touching:
            if not other.classification == vertex.classification:
                rim = True
                continue         
            for patch in patches:
                if patch.classification == vertex.classification and other in patch.vertices:
                    patchable.add(patch)
        new_patch = Patch(vertex, rim)
        if vertex.highest:
            new_patch.highest = True
        patches.add(new_patch)
        for patch in patchable:
            new_patch.add_patch(patch)
            patches.remove(patch)
    checked = set()
    d = 0
    for patch in patches:
        if not patch.core:
            d += 1 
            continue
        patch.finish()
        checked.add(patch)
    return checked, d

def Build_Patch_Distances(distances_dict: dict, patches: Set[Patch]) -> None:
        for patch in patches:
            distances = patch.distances
            for rim in patch.rim:
                number = rim.number
                touching = sortedset()
                for vertex in rim.touching:
                    if vertex.classification == rim.classification:
                        touching.append(vertex, distances_dict[frozenset((rim, vertex))])
                while touching:
                    vertex, value = touching.pop(0)
                    distances[vertex.number][number] = value
                    for other in vertex.touching:
                        if not other.classification == rim.classification:
                            if other.patch:
                                other.patch.touching.add(rim.patch)
                                rim.patch.touching.add(other.patch)
                            continue
                        if other is rim or patch.distances[other.number][number]:
                            continue
                        touching.append(other, value + distances_dict[frozenset((vertex, other))])

def Build_Surfaces(patches: Set[Patch]) -> Tuple[Set[FrozenSet[Patch]], FrozenSet[Patch]]:
    unassagined = patches.copy()
    inner = set()
    outer = None
    while unassagined:
        patch = unassagined.pop()
        surface = {patch}
        touching = patch.touching.copy()
        highest = patch.highest
        while touching:
            patch = touching.pop()
            surface.add(patch)
            unassagined.remove(patch)
            touching |= patch.touching & unassagined
            if patch.highest: highest = True
        if highest:
            outer = frozenset(surface)
        else:
            inner.add(frozenset(surface))
    for surface in inner:
        for patch in surface:
            patch.surface = surface
    for patch in outer:
        patch.surface = outer
    return inner, outer

def Get_Centers(patches: Set[Patch]) -> None:
    for patch in patches:
        patch.get_center()

def Build_Distances_Dict(distances_dict: dict, patches: Set[Patch]) -> Set[FrozenSet]:
    unchecked = patches.copy()
    patch = unchecked.pop()
    while unchecked:
        to_check = unchecked & patch.surface
        if not to_check:
            patch = unchecked.pop()
            continue
        touching = sortedset()
        if patch.center:
            center = patch.center
            for vertex in center.touching:
                touching.append(vertex, distances_dict[frozenset((center, vertex))])
        while touching:
            vertex = touching.pop(0)[0]
            for other in vertex.touching:
                if other is center:
                    continue
                center_to_other = frozenset((center, other))
                if center_to_other in distances_dict:
                    if other.center:
                        touching.append(other, distances_dict[center_to_other])
                    continue
                distance = distances_dict[frozenset((center, vertex))] + distances_dict[frozenset((vertex, other))]
                distances_dict[center_to_other] = distance
                touching.append(other, distance)
                if not other.center:
                    continue
                to_check.remove(other.patch)
                if not to_check:
                    patch = other.patch
                    unchecked.remove(patch)
                    break
            else:
                continue
            break

def Surface_Analysis(protein: str, logger, start) -> None:
    file = protein+"/"+protein+'.wrl'
    distances_dict = dict()
    end = time()
    logger.info('Started analizing protein.'.format(protein), extra={'offset': timedelta(seconds=end-start), 'code': protein})
    vertices = Read_Wrl(distances_dict, file)
    end = time()
    logger.debug('Protein model imported from pymol.', extra={'offset': timedelta(seconds=end-start), 'code': protein})
    patches, d = Get_Patches(vertices)
    end = time()
    logger.debug(str(len(patches))+' patches built, '+str(d)+' deleted due to small size.', extra={'offset': timedelta(seconds=end-start), 'code': protein})
    Build_Patch_Distances(distances_dict, patches)
    end = time()
    logger.debug('Internal patch distances calculated.', extra={'offset': timedelta(seconds=end-start), 'code': protein})
    inners, outer = Build_Surfaces(patches)
    end = time()
    logger.debug(str(len(inners)+1)+' surfaces built.', extra={'offset': timedelta(seconds=end-start), 'code': protein})
    Get_Centers(patches)
    end = time()
    logger.debug('Patch centers found.', extra={'offset': timedelta(seconds=end-start), 'code': protein})
    Build_Distances_Dict(distances_dict, patches)
    end = time()
    logger.debug('Distances between patches calculated.', extra={'offset': timedelta(seconds=end-start), 'code': protein})
    try:
        sizes_outer_file = open(protein+'/sizes_outer_'+protein+'.csv'.format(protein), 'w', newline='')
        sizes_outer_writer = writer(sizes_outer_file)
        sizes_outer_writer.writerow(['Classification', 'Size(\u212B\u00B2)'])
        sizes_inner_file = open(protein+'/sizes_inner_'+protein+'.csv'.format(protein), 'w', newline='')
        sizes_inner_writer = writer(sizes_inner_file)
        sizes_inner_writer.writerow(['Classification', 'Size(\u212B\u00B2)', 'Surface'])
        distances_outer_file = open(protein+'/distances_outer_'+protein+'.csv', 'w', newline='')
        distances_outer_writer = writer(distances_outer_file)
        distances_outer_writer.writerow(['Classification', 'Distance(\u212B)'])
        distances_inner_file = open(protein+'/distances_inner_'+protein+'.csv', 'w', newline='')
        distances_inner_writer = writer(distances_inner_file)
        distances_inner_writer.writerow(['Classification', 'Distance(\u212B)', 'Surface'])
        touching_outer_file = open(protein+'/touching_outer_'+protein+'.csv', 'w', newline='')
        touching_outer_writer = writer(touching_outer_file)
        touching_outer_writer.writerow(['First Classification', 'Second Classification', 'Distance(\u212B)'])
        touching_inner_file = open(protein+'/touching_inner_'+protein+'.csv', 'w', newline='')
        touching_inner_writer = writer(touching_inner_file)
        touching_inner_writer.writerow(['First Classification', 'Second Classification', 'Distance(\u212B)', 'Surface'])
        for i, surface in enumerate(inners):
            checked = set()
            for first in surface:
                sizes_inner_writer.writerow([first.classification, first.area(), i])
                for second in surface:
                    if second is first or second in checked:
                        continue
                    if second.classification == first.classification:
                        distances_inner_writer.writerow([first.classification, distances_dict[frozenset((first.center, second.center))], i])
                    if first.is_touching_patch(second):
                        touching_inner_writer.writerow([first.classification, second.classification, distances_dict[frozenset((first.center, second.center))], i])
                    checked.add(first)
        checked = set()
        for first in outer:
            sizes_outer_writer.writerow([first.classification, first.area()])
            for second in outer:
                if second is first or second in checked:
                    continue
                if second.classification == first.classification:
                    distances_outer_writer.writerow([first.classification, distances_dict[frozenset((first.center, second.center))]])
                if first.is_touching_patch(second): 
                    touching_outer_writer.writerow([first.classification, second.classification, distances_dict[frozenset((first.center, second.center))]])
            checked.add(first)
        open(protein+'/.finished', 'w')
    except Exception as e:
        raise e
    finally:
        sizes_outer_file.close()
        sizes_inner_file.close()
        distances_outer_file.close()
        distances_inner_file.close()
        touching_outer_file.close() 
        touching_inner_file.close()   
    end = time()
    logger.info('Finished analizing protein.', extra={'offset': timedelta(seconds=end-start), 'code': protein})