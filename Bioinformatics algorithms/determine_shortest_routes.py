def bepaalAantalKortsteRoutes(n, m, openP):
    road = {(0,0):1}
    points = 0
    for row in range(0,n):
        for col in range(0,m):
            if (row,col) in openP:
                road[(row, col)] = 0
            if (row,col) not in road:
                if col != 0:
                    points += road[(row, col-1)]
                if row != 0:
                    points += road[(row-1, col)]
                road[(row, col)] = points
                points = 0
    return road[n-1, m-1]