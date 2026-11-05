from objects.vector import Vector
from objects.line import Line
from objects.point import Point
import math


class Ray:
    """
    Class representing the actual sound Ray in the room. Reflects of the walls with correct angles.
    The order given at the start is the amount of times this specific ray is allowed to be reflected.
    """
    def __init__(self, r: Vector, order: int):
        assert order >= 0
        self.value = r
        self.values = []
        self.expanded = False
        self.order = order
        self.recieved = False # Flag to mark, if we recieved the ray, meaning we dont want to continue anymore.
        self.logger = []

    # ------------------------------------------------------------------------

    def expand(self, walls: list[Line], reciever: tuple[Point, float], show: bool = False) -> list[Vector]:
        """
        Expands the given vector, until it was reflected order times.
        
        :param self: Starting vector
        :param walls: The walls, where the ray gets reflected off.
        :type walls: list[Line]
        :param reciever: The position of the reciever and the radius, which we count as recieved in.
        :type reciever: tuple[Point, float]
        :param show: Add all the happening stuff to the plot (intersection, imige method, etc.)
        :type show: bool
        :return: Returns and saves a list of the vectors, which are reflection from the starting point.
        :rtype: list[Vector]
        """
    
        # If this ray already got expanded, just return the saved one, because environment is static.
        if self.expanded:
            return self.values
        self.expanded = True

        prev = None # saving the previous intersection point, because it will be the new anchor.
        block = -1 # flag for blocking reflection of the same wall twice in a row.
        self.value = self.value.extend() # we want a line, to perform point intersection in the first iteration.
        
       

        # repeating the process for each order of reflection.
        for o in range(0, self.order + 1):
            curr = self.value

            # saving the intersection with each wall.
            intersections = []
            for j, w in enumerate(walls):
                if j != block:

                    # compute the intersection with the wall.
                    intersections.append(w.intersection(curr))
                else:
                    intersections.append(None)
            
            # intersection with reciever.
            intersections.append(curr.intersection(reciever))

            # if recieved, this is the last iteration.
            if not intersections[-1] is None:
                self.recieved = True
                
            
            # saving each intersecting vector and its distance as a deciding metric.
            vecs = []
            dist = []
            for i, inter in enumerate(intersections):

                # If there was an intersection
                if not inter is None:

                    # create the vector measure its distance.
                    vecs.append(Vector(curr.anchor, inter))
                    dist.append(abs(vecs[i]))

                    # guarding a special case, because lines/walls are infinit and a reflection can happen with
                    # a wall ouside of the body.
                    # There for the distance needs to be atleast to the previous intersection, to enter the 
                    # figure again.
                    if (not prev is None) and abs(vecs[i]) <= abs(Vector(curr.anchor, prev)):
                        dist[i] = math.inf
                
                # No interscection, then add infinity, so it is overlooked.
                else:
                    vecs.append(None)
                    dist.append(math.inf)
            
            # take the closest legal intersecton
            min_d = min(dist)
            if min_d == math.inf:
                break

            m = dist.index(min_d)
            block = m
            self.logger.append(m + 1)
            
            # adjust the respective line.
            if self.recieved:
                int_point = intersections[-1]
                new_line = None
            
            else:
                int_point = walls[m].intersection(curr)
                if int_point is None:
                    break
                new_line = walls[m].image(int_point, curr.anchor, show, o) # the line
        

            # for the first case
            if prev is None:
                self.values.append(Vector(curr.anchor, int_point))
            
            # for all the other cases.
            else:
                self.values.append(Vector(prev, int_point))
            

            prev = int_point
            self.value = new_line

            # we recieved the signal in this iteration, so we must not continue.            
            if self.recieved:
                break
        
        # if we are done, we return and save our computed values.
        self.expanded = True
        return self.values
    
    # ------------------------------------------------------------------------

    def plot(self, color: str = "black", name: str = ""):
        """
        Plots the ray with all its reflections, if they exist.
        Uses the vectors plotting method.
        
        :param self: Collected values.
        :param color: color, default is black
        :type color: str
        :param name: label for the ray to be displayed.
        :type name: str
        """
        for value in self.values:
            value.plot(color, name)
    
    # ------------------------------------------------------------------------

    def pov(self, reciever: tuple[Point, float], color: str = "orange", name: str = "", labelling: bool = True) -> tuple[Vector, str]:
        """
        Plotting every Point, where the reciever sees it, just using the incoming direction
        and the entire length of the ray for the plot.
        Only works, if this specific ray object is recieved.
        
        :param self: All the previous and the last reflection vectors.
        :param reciever: The coordinates of the reciever and the radius as a buffer.
        :param color: color for the plot, with default orange
        :type color: str
        :param name: label
        :type name: str
        """
        rec, rad = reciever

        # only continue if this ray was ever recieved        
        if not self.recieved:
            return None
        
        # entire distance of the ray
        distance = sum(abs(value) for value in self.values)

        # direction of intersection
        pov_vec = Vector(rec, self.values[-1].anchor)

        # normalized
        npov = pov_vec.normalized()

        # scaled
        pov = npov.scalar(distance + rad)

        # adding the label:
        
        if labelling:
            label = "I"
            for reflection in self.logger:
                label += str(reflection)
            label = label[:len(label) - 1]
        
        pov.end.plot(color, label)
        
        
        return (pov, label)


