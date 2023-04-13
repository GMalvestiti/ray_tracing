# -*- coding: utf-8 -*-

class Screen(object):

    # Construtor da classe
    def __init__(self, title, bgColor, width, height):
        self.title = title       # título da janela
        self.bgColor = bgColor   # cor de fundo
        self.width = width       # largura da janela
        self.height = height     # altura da janela
        self.screen = pygame.display.set_mode(self.size()) # define o tamanho da tela
        pygame.display.set_caption(self.title)             # define o título da janela
        self.clock = pygame.time.Clock()
        
        
    # Executa o pipeline gráfico
    def run(self, obj):
        while True:  # laço principal
            # captura eventos
            for event in pygame.event.get(): 
                
                # Captura evento de clicar em botão para fechar
                if event.type == pygame.QUIT:
                    return pygame.quit()
            
            # preencha a tela com a cor de fundo
            self.screen.fill(self.bgColor)
            
            # gera o desenho
            obj.draw(self)
            
            # aplica o antialiasing
            # self.meanFilter()
            
            # atualiza a tela 
            pygame.display.update()
            
            self.clock.tick(30)
        
        
    # retorna um vetor com o tamanho da tela
    def size(self):
        return (self.width, self.height)
    
    # modifica um pixel na tela com a cor desejada
    def setPixel(self, x, y, color):
        self.screen.set_at((x, y), color)
    
    # filtro da média para o antialising
    def meanFilter(self):
        # Captura a matrix da tela
        #frameBuffer2 = pygame.PixelArray(self.screen)
        
        from copy import copy
        frameBuffer = pygame.surfarray.array3d(self.screen)
        #print(frameBuffer)
        
        import numpy as np
        mask = np.ones((3, 3)) * 1/9 
        
        #print(mask)
               
        for i in range(1, self.width - 1):
            for j in range(1, self.height - 1):               
                temp = np.zeros((3))
                
                for k in range(-1,2):
                    for l in range(-1,2):
                        for b in range(3):
                            temp[b] = temp[b] + frameBuffer[i + k][j + l][b] * mask[k + 1][l + 1]
                        
                #print(pygame.Color(int(temp[0]), int(temp[1]), int(temp[2]), 255))
                        
                self.setPixel(i, j, pygame.Color(int(temp[0]), int(temp[1]), int(temp[2]), 255));
                #frameBuffer2[i][j] = pygame.Color(int(temp[0]), int(temp[1]), int(temp[2]), 255);

class CoordHomog(object):
    def __init__(self, x, y, z, w):
        self.x = x
        self.y = y
        self.z = z
        self.w = w
        
    @property
    def x(self):
        return self._x
    
    @x.setter
    def x(self, x):
        self._x = x
    
    @property
    def y(self):
        return self._y
    
    @y.setter
    def y(self, y):
        self._y = y
        
    @property
    def z(self):
        return self._z
    
    @z.setter
    def z(self, z):
        self._z = z
    
    @property
    def w(self):
        return self._w
    
    @w.setter
    def w(self, w):
        self._w = w
        
    @staticmethod
    def getCoordFromArray(v):
        return CoordHomog(v[0], v[1], v[2], v[3])
        
    def getArray(self):
        return np.array([self.x, self.y, self.z, self.w])
    
    def getRealCoord(self):
        return np.array([self.x / self.w, self.y / self.w, self.z / self.w, 1])
    
    def divideByHomogeneousCoord(self):
        self.x = self.x / self.w
        self.y = self.y / self.w
        self.z = self.z / self.w
        self.w = 1.0

# Definição da Classe para Primitiva de Linha

import numpy as np

class Primitive(object):
    # construtor da classe
    def __init__(self, id, type):
        self.id = id
        self.type = type

class Line(Primitive):
    # construtor da classe
    def __init__(self, p1, p2, color, id):
        super().__init__(id, 'line')
        self.p1 = p1       # coordenada x do primeiro ponto
        self.p2 = p2       # coordenada x do segundo ponto
        self.color = color # cor do objeto
    
    def getPoint(self):
        return np.array([self.p1.getArray(), self.p2.getArray()])
    
    # renderiza a linha desejada na tela
    def draw(self, screen):
        self.dda(screen)

    # Algoritmo DDA
    def dda(self, screen):      
        # Definição e Inicialização de Variáveis locais
        dx, dy, k = 0, 0, 0
        x_inc, y_inc = 0.0, 0.0
        x, y = 0.0, 0.0
        x2 = int(round(self.p2.x))
        y2 = int(round(self.p2.y))
        x1 = int(round(self.p1.x))
        y1 = int(round(self.p1.y))
    
        # Define os deslocamentos nas direções x e y
        dx = x2 - x1
        dy = y2 - y1
    
        # Define qual a direção de incremento fixo
        if abs(dx) > abs(dy):
            iter = abs(dx)
        else:
            iter = abs(dy)
            
        # Define o ponto inicial
        x = x1
        y = y1
        
        # Desenha o ponto inicial na tela
        screen.setPixel(round(x), round(y), self.color)
        
        if iter == 0:
            return
        
        # Define os incrementos para cada direção
        x_inc = dx/iter
        y_inc = dy/iter

        # Geração e renderização dos pontos seguintes da linha
        for k in range(iter):
            # Gera o próximo ponto
            x = x + x_inc
            y = y + y_inc
            
            # Desenha o ponto
            screen.setPixel(round(x), round(y), self.color)
            
    def bresenham(self, screen):
        # Definição e inicialização de variáveis locais
        dx, dy, d = 0, 0, 0
        incrE, incrNE = 0, 0
        x, y, xFinal = 0, 0, 0
        
        # Define os deslocamentos absolutos nas direções x e y
        dx = abs(self.p2.x - self.p1.x)
        dy = abs(self.p2.y - self.p1.y)
        
        # Define o d de teste inicial
        d = 2 * dy - dx
        
        # Define os incrementos nas direções x e y
        incrE = 2 * dy
        incrNE = 2 * (dy - dx)
        
        # Troca a ordem dos pontos em caso de segundo ponto à esquerda de primeiro ponto
        if self.p1.x > self.p2.x:
            x = self.p2.x
            y = self.p2.y
            xFinal = self.p1.x
        else:
            x = self.p1.x
            y = self.p1.y
            xFinal = self.p2.x
        
        # Desenha o ponto inicial na tela
        screen.setPixel(x, y, self.color)
        
        # Gera e renderiza os pontos seguintes da linha
        while x < xFinal:
            # Gera o próximo ponto
            x = x + 1
            
            if d < 0:
                d = d + incrE
            else:
                y = y + 1
                d = d + incrNE
            
            # Desenha o próximo ponto
            screen.setPixel(x, y, self.color)
    
    def clip(self, window):
        return self.cohen_sutherland(window)
    
    #def clip3D(self, window, near, far):
    #    return self.cohen_sutherland3D(window, near, far)
    
    def cohen_sutherland(self, window):
        xI = self.p1.x
        yI = self.p1.y
        xF = self.p2.x
        yF = self.p2.y
        
        # Cria o código de região dos pontos
        p1Code = self.codify(xI, yI, window)
        p2Code = self.codify(xF, yF, window)
        
        while True:
            # Testa inicial se aceita ou rejeita
            # Totalmente dentro?
            if p1Code | p2Code == 0: # totalmente dentro -> ACEITA
                pI = Point3D(xI, yI, 0, 1, self.p1.color, self.p1.id)
                pF = Point3D(xF, yF, 0, 1, self.p2.color, self.p2.id)
                return [Line(pI, pF, self.color, self.id)]
        
            # Totalmente fora?
            if p1Code & p2Code != 0: # totalmente fora -> DESCARTA
                return []
        
            # Há intersecção
            # Verifica contra todos os limites
            if p1Code != 0:
                pTest = p1Code
            else:
                pTest = p2Code
            
            if pTest & 1: # ESQ      
                x = window.xMin
                y = yI + (yF - yI) / (xF - xI) * (window.xMin - xI)
            elif pTest & 2: # DIR
                x = window.xMax
                y = yI + (yF - yI) / (xF - xI) * (window.xMax - xI)
            elif pTest & 4: # ABAIXO
                x = xI + (window.yMin - yI) * (xF - xI) / (yF - yI)
                y = window.yMin
            elif pTest & 8: # ACIMA
                x = xI + (window.yMax - yI) * (xF - xI) / (yF - yI)
                y = window.yMax
            
            if pTest == p1Code:
                xI = x
                yI = y
                p1Code = self.codify(xI, yI, window)
            else:
                xF = x
                yF = y
                p2Code = self.codify(xF, yF, window)
        
    def codify(self, x, y, window):
        return np.signbit(window.yMax - y) * 8 + np.signbit(y - window.yMin) * 4 + np.signbit(window.xMax - x) * 2 + np.signbit(x - window.xMin)
    
    def transform(self, T):
        point = self.p1.getArray().T
        res = T.matrix.dot(point)
        #print(res)
        
        point2 = self.p2.getArray().T
        res2 = T.matrix.dot(point2)
        
        return Line(Point3D(res[0], res[1], res[2], res[3], self.p1.color, self.p1.id), Point3D(res2[0], res2[1], res2[2], res2[3], self.p2.color, self.p2.id), self.color, self.id)
    
    def genRealPoints(self):
        self.p1.divideByHomogeneousCoord()
        self.p2.divideByHomogeneousCoord()

# Programa Principal

# Carregamento de bibliotecas
import pygame

# Inicialização do PyGame
pygame.init()

# Definição da Classe para Primitiva de Circunferência

class Circle(Primitive):
    # construtor da classe
    def __init__(self, xc, yc, raio, color, id):
        super().__init__(id, 'circle')
        self.xc = xc       # coordenada x do centro
        self.yc = yc       # coordenada y do centro
        self.raio = raio   # raio da circunferência
        self.color = color # cor do objeto
    
    # renderiza a circunferência desejada na tela
    # Algoritmo de Bresenham
    def draw(self, screen):
        # Definição e Inicialização de Variáveis locais      
        # Coniderando a circunferência ao redor da origem, mas renderizada transladada
        x = 0                # x inicial
        y = self.raio        # y inicial
        d = 1 - self.raio    # d de teste inicial

        # Desenha os pontos inicias de cada quadrante
        self.drawCirclePoints(x, y, screen)
        
        # Gera os novos pontos e os renderiza
        while x < y:
            
            if d < 0:   # Direção E
                d = d + 2 * x + 3
            else:       # Direção SE
                d = d + 2 * (x - y) + 5
                y = y - 1
            
            x = x + 1

            self.drawCirclePoints(x, y, screen)
            
    def drawCirclePoints(self, x, y, screen):
        xCentro = self.xc
        yCentro = self.yc
        screen.setPixel(xCentro + x, yCentro + y, self.color)
        screen.setPixel(xCentro + y, yCentro + x, self.color)
        screen.setPixel(xCentro + y, yCentro - x, self.color)
        screen.setPixel(xCentro + x, yCentro - y, self.color)
        screen.setPixel(xCentro - x, yCentro - y, self.color)
        screen.setPixel(xCentro - y, yCentro - x, self.color)
        screen.setPixel(xCentro - y, yCentro + x, self.color)
        screen.setPixel(xCentro - x, yCentro + y, self.color)

# Criação de uma classe que contém o desenho

class Picture(object):
    # Construtor da Classe
    def __init__(self):
        self.primitivas = [] # Define uma lista de primitivas para representar um desenho
    
    def draw(self, screen):
        # Telhado
        p1 = CoordHomog(350,50, 0, 1)
        p2 = CoordHomog(50, 250, 0, 1)
        p3 = CoordHomog(650, 250, 0, 1)
        l1 = Line(p1, p2, pygame.Color(255, 0, 0, 255), 1)
        l2 = Line(p2, p3, pygame.Color(255, 0, 0, 255), 2)
        l3 = Line(p3, p1, pygame.Color(255, 0, 0, 255), 3)
        
        # Parede
        p4 = CoordHomog(50, 650, 0, 1)
        p5 = CoordHomog(650, 650, 0, 1)
        l4 = Line(p2, p4, pygame.Color(0, 200, 100, 255), 4)
        l5 = Line(p4, p5, pygame.Color(0, 200, 100, 255), 5)
        l6 = Line(p5, p3, pygame.Color(0, 200, 100, 255), 6)
        
        circ = Circle(350, 450, 150, pygame.Color(0, 0, 255, 255), 7);
        
        # Insere primitivas na lista
        self.primitivas.append(l1);
        self.primitivas.append(l2);
        self.primitivas.append(l3);
        self.primitivas.append(l4);
        self.primitivas.append(l5);
        self.primitivas.append(l6);
        self.primitivas.append(circ);
        
        # Desenha cada primitiva que está na lista
        for item in self.primitivas:
            item.draw(screen);

# Classe da primitiva Ponto

class Point2D(object):
    # construtor da Classe
    def __init__(self, x, y, color):
        self.x = x
        self.y = y
        self.color = color
    
    # Renderiza um ponto
    def draw(self, screen):
        screen.setPixel(self.x, self.y, self.color)
        
    @property
    def x(self):
        return self._x
    
    @x.setter
    def x(self, x):
        self._x = x
    
    @property
    def y(self):
        return self._y
    
    @y.setter
    def y(self, y):
        self._y = y
        
    @property
    def color(self):
        return self._color
    
    @color.setter
    def color(self, color):
        self._color = color

class Point3D(Primitive, CoordHomog):
    # construtor da Classe
    def __init__(self, x, y, z, w, color, id):
        Primitive.__init__(self, id, 'point')
        CoordHomog.__init__(self, x, y, z, w)
        self.color = color
            
    @property
    def color(self):
        return self._color
    
    @color.setter
    def color(self, color):
        self._color = color
        
    # Renderiza um ponto
    def draw(self, screen):
        screen.setPixel(int(round(self.x)), int(round(self.y)), self.color)

# Classe para informação de Arestas

class EdgeInfo(object):
    # Construtor da Classe
    def __init__(self, initialPoint, finalPoint):
        if int(round(initialPoint.y)) <= int(round(finalPoint.y)):
            self.yMax = int(round(finalPoint.y))
            self.x = int(round(initialPoint.x))  # x corrente, inicialmente x in Ymin
            self.yMin = int(round(initialPoint.y))
        else:
            self.yMax = int(round(initialPoint.y))
            self.x = int(round(finalPoint.x))     # x corrente, inicialmente x in Ymin
            self.yMin = int(round(finalPoint.y))
            
        self.inverseOfAngularCoefficient = float(round(finalPoint.x) - round(initialPoint.x)) \
                                            / float(round(finalPoint.y) - round(initialPoint.y))
            
    @property
    def yMax(self):
        return self._yMax
    
    @yMax.setter
    def yMax(self, yMax):
        self._yMax = yMax
    
    @property
    def yMin(self):
        return self._yMin
    
    @yMin.setter
    def yMin(self, yMin):
        self._yMin = yMin
    
    @property
    def x(self):
        return self._x
    
    @x.setter
    def x(self, x):
        self._x = x
    
    def updateX(self):
        self.x = self.x + self.inverseOfAngularCoefficient

# Classe da primitiva Polígono

class Polygon(Primitive):
    # construtor da classe
    def __init__(self, showEdges, edgeColor, isFilled, fillColor, id):
        super().__init__(id, 'polygon')
        self.listOfPoints = []
        self.showEdges = showEdges
        self.edgeColor = edgeColor
        self.isFilled = isFilled
        self.fillColor = fillColor
    
    # Adiciona vértices na lista
    def addVertex(self, point):
        self.listOfPoints.append(point)
    
    # Renderiza Polígono
    def draw(self, screen):
        if len(self.listOfPoints) < 3:
            print("Não forma polígono. Menos de 3 vértices.")
            return
            
        # Desenha arestas se desejar
        if self.isFilled:
            self.scanline(screen)
        
        if self.showEdges:
            for i in range(0,len(self.listOfPoints) - 1):     
                pI = self.listOfPoints[i]
                pF = self.listOfPoints[i+1]
                line = Line(pI, pF, self.edgeColor, i)
                line.draw(screen)
            
            pI = self.listOfPoints[-1]
            pF = self.listOfPoints[0]
            line = Line(pI, pF, self.edgeColor, len(self.listOfPoints) - 1)
            line.draw(screen)
            
            for i in range(0,len(self.listOfPoints)):               
                self.listOfPoints[i].draw(screen)
    
    # Faz o scanline para preencher o polígono
    def scanline(self, screen): 
        z_buffer = []
        frameBuffer = []
        width = screen.size()[0]
        height = screen.size()[1]
        for i in range(height):
            linha = []
            linha2 = []
            for j in range(width):
                linha.append(1.0)
                linha2.append(screen.color())
            z_buffer.append(linha)
            frameBuffer.append(linha2)
        
        yMax = int(round(self.listOfPoints[0].y))
        for item in self.listOfPoints:
            if int(round(item.y)) > yMax:
                yMax = int(round(item.y))
                
        x2 = item.x() - item.inverseOfAngularCoefficient()
                
        y = yMax # armazena o y corrente, começando pelo valor mínimo
        
        #### Cria tabela de arestas ####
        edgeTable = []
        for i in range (0, yMax+1):
            edgeTable.append([])
            
        for i in range(0,len(self.listOfPoints) - 1):     
        #for i in range(0,1):     
            # exclui arestas horizontais
            if round(self.listOfPoints[i].y) - round(self.listOfPoints[i+1].y) != 0:
                edge = EdgeInfo(self.listOfPoints[i], self.listOfPoints[i+1])
                yMin = edge.yMin
                if yMin < y:
                    y = yMin
    
                print(yMin)
                edgeTable[yMin].append(edge)       
        
        # Fecha o polígono
        # exclui arestas horizontais
        if round(self.listOfPoints[-1].y) - round(self.listOfPoints[0].y) != 0:
            edge = EdgeInfo(self.listOfPoints[-1], self.listOfPoints[0])
            yMin = edge.yMin
            if yMin < y:
                y = yMin
            edgeTable[yMin].append(edge)
      
        ####
        activeET = []
        
        ### Laço principal
        while y <= yMax:
            lp = self.getListOfPoints();
            z = ((-1)*(lp[0])*(x2+1)-lp[0]*y-lp[3])/lp[2]
            z2 = z - lp[0]/lp[2]
            if z < z_buffer[x][y]:
                z_buffer[x][y] = z
                frameBuffer[x][y] = screen.color()
                    
            # Move a lista y na ET para AET (ymin = y), mantendo a AET ordenada em x
            activeET.extend(edgeTable[y])
            edgeTable[y] = []
            activeET.sort(key = sortByX)
            
            # Desenhe os pixels do bloco na linha de varredura y, 
            # usando os pares de coordenadas x da AET (cada dois nós definem um bloco)
            for i in range(0, len(activeET) - 1, 2):
                for x in range(int(activeET[i].x), int(activeET[i + 1].x + 1)):
                    screen.setPixel(x, y, self.fillColor)
            
            # Atualiza o valor de y para a próxima linha de varredura
            y = y + 1
            
            # Remova as arestas que possuem ymax = y da AET
            delL = []
            for item in activeET:
                if item.yMax <= y:
                    delL.append(item)
            
            for item in delL:
                activeET.remove(item)
            
            delL.clear()
            
            # Para cada aresta na AET, atualize x = x + 1/m
            for item in activeET:
                item.updateX()
                
                
    def empty(self, ET):
        for item in ET:
            if item:
                return False
        
        return True
    
    def getListOfPoints(self):
        return self.listOfPoints
    
    def transform(self, T):
        p = Polygon(self.showEdges, self.edgeColor, self.isFilled, self.fillColor, self.id)
        for item in self.listOfPoints:
            point = item.getArray().T
            #print(T.matrix)
            res = T.matrix.dot(point)
            #print(res)
            p.addVertex(Point3D(res[0], res[1], res[2], res[3], item.color, item.id))
        
        return p
    
    def genRealPoints(self):
        for item in self.listOfPoints:
            item.divideByHomogeneousCoord()
        
    def sutherland_hodgman(self, window):
        p = self              
        
        for k in range (1,5):
            temp = Polygon(self.showEdges, self.edgeColor, self.isFilled, self.fillColor, self.id)
            #print(len(p.listOfPoints))
            for i in range(1,len(p.listOfPoints) + 1):     
                #print(i - 1)
                #print(len(p.listOfPoints))
                pI = p.listOfPoints[i - 1]
                pF = p.listOfPoints[i % len(p.listOfPoints)]
                
                if k == 1: # ESQ
                    if (pI.x >= window.xMin):
                        if (pF.x >= window.xMin): # Ambos dentro
                            temp.addVertex(pF)
                        else: # Dentro para Fora
                            m = (pF.y - pI.y) / (pF.x - pI.x)
                            yIntersected = pI.y + m * (window.xMin - pI.x)
                            temp.addVertex(Point3D(window.xMin, yIntersected, 0, 1, pF.color, pF.id))
                    else:
                        if (pF.x >= window.xMin): # Fora para Dentro
                            m = (pF.y - pI.y) / (pF.x - pI.x)
                            yIntersected = pI.y + m * (window.xMin - pI.x)
                            temp.addVertex(Point3D(window.xMin, yIntersected, 0, 1, pF.color, pF.id))
                            temp.addVertex(pF)
                            
                elif k == 2: # DIR
                    if (pI.x <= window.xMax):
                        if (pF.x <= window.xMax): # Ambos dentro
                            temp.addVertex(pF)
                        else: # Dentro para Fora
                            m = (pF.y - pI.y) / (pF.x - pI.x)
                            yIntersected = pI.y + m * (window.xMax - pI.x)
                            temp.addVertex(Point3D(window.xMax, yIntersected, 0, 1, pF.color, pF.id))
                    else:
                        if (pF.x <= window.xMax): # Fora para Dentro
                            m = (pF.y - pI.y) / (pF.x - pI.x)
                            yIntersected = pI.y + m * (window.xMax - pI.x)
                            temp.addVertex(Point3D(window.xMax, yIntersected, 0, 1, pF.color, pF.id))
                            temp.addVertex(pF)
                            
                elif k == 3: # ABAIXO
                    if (pI.y >= window.yMin):
                        if (pF.y >= window.yMin): # Ambos dentro
                            temp.addVertex(pF)
                        else: # Dentro para Fora
                            mInv = (pF.x - pI.x) / (pF.y - pI.y)
                            xIntersected = pI.x + (window.yMin - pI.y) * mInv
                            temp.addVertex(Point3D(xIntersected, window.yMin, 0, 1, pF.color, pF.id))
                    else:
                        if (pF.y >= window.yMin): # Fora para Dentro
                            mInv = (pF.x - pI.x) / (pF.y - pI.y)
                            xIntersected = pI.x + (window.yMin - pI.y) * mInv
                            temp.addVertex(Point3D(xIntersected, window.yMin, 0, 1, pF.color, pF.id))
                            temp.addVertex(pF)
                elif k == 4: # ACIMA
                    if (pI.y <= window.yMax):
                        if (pF.y <= window.yMax): # Ambos dentro
                            temp.addVertex(pF)
                        else: # Dentro para Fora
                            mInv = (pF.x - pI.x) / (pF.y - pI.y)
                            xIntersected = pI.x + (window.yMax - pI.y) * mInv
                            temp.addVertex(Point3D(xIntersected, window.yMax, 0, 1, pF.color, pF.id))
                    else:
                        if (pF.y <= window.yMax): # Fora para Dentro
                            mInv = (pF.x - pI.x) / (pF.y - pI.y)
                            xIntersected = pI.x + (window.yMax - pI.y) * mInv
                            temp.addVertex(Point3D(xIntersected, window.yMax, 0, 1, pF.color, pF.id))
                            temp.addVertex(pF)
            p = temp
                
        return [p]
    
    def clip(self, window):
        if self.isFilled:
            return self.sutherland_hodgman(window)
        else:
            p = []
            for i in range(0,len(self.listOfPoints) - 1):     
                pI = self.listOfPoints[i]
                pF = self.listOfPoints[i+1]
                line = Line(pI, pF, self.edgeColor, i)
                p.extend(line.clip())
            
            pI = self.listOfPoints[-1]
            pF = self.listOfPoints[0]
            line = Line(pI, pF, self.edgeColor, len(self.listOfPoints) - 1)
            p.extend(line.clip())
            
            return p
                

# Usada para ordenar a AET por valores de x
def sortByX(item):
    return item.x

import numpy as np
import math

class GeometricTransformation(object):
    # construtor da Classe
    def __init__(self, matrix):
        self.matrix = matrix

    @property
    def matrix(self):
        return self._matrix
    
    @matrix.setter
    def matrix(self, matrix):
        self._matrix = matrix


class Translation3D(GeometricTransformation):
    # construtor da classe
    def __init__(self, Tx, Ty, Tz):
        super().__init__(np.array([[1, 0, 0, Tx], [0, 1, 0, Ty], [0, 0, 1, Tz], [0, 0, 0, 1]]))
    
class Scale3D(GeometricTransformation):
    # construtor da classe
    def __init__(self, Sx, Sy, Sz):
        super().__init__(np.array([[Sx, 0, 0, 0], [0, Sy, 0, 0], [0, 0, Sz, 0], [0, 0, 0, 1]]))
        
class Rotation3D(GeometricTransformation):
    # construtor da classe
    def __init__(self, theta, axis):
        self.axis = axis
        
        if axis == 'z':
            super().__init__(np.array([[math.cos(theta), -math.sin(theta), 0, 0], [math.sin(theta), math.cos(theta), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]))
        elif axis == 'x':
            super().__init__(np.array([[1, 0, 0, 0], [0, math.cos(theta), -math.sin(theta), 0], [0, math.sin(theta), math.cos(theta), 0], [0, 0, 0, 1]]))
        elif axis == 'y':
            super().__init__(np.array([[math.cos(theta), 0, math.sin(theta), 0], [0, 1, 0, 0], [-math.sin(theta), 0, math.cos(theta), 0], [0, 0, 0, 1]]))
        else:
            print("Eixo não definido!")
            exit()
            
        #print(self.matrix)

class CombinedTransformation(GeometricTransformation):
    # construtor da classe
    def __init__(self):
        super().__init__([])
        self.stack = []
    
    def add(self, T):
        self.stack.append(T)
    
    #def pop(self):
    #    return self.stack.pop(-1)
              
    def combine(self):
        while (len(self.stack) > 1):
            m2 = self.stack.pop(-1)
            #print(m2.matrix)
            m1 = self.stack.pop(-1)
            #print(m1.matrix)
            m3 = GeometricTransformation(m2.matrix.dot(m1.matrix))
            #print(m3.matrix)
            
            #print('Novo')
            self.stack.append(m3)
            
        # Retira o item restante da pilha    
        # Esta é a matrix combinada
        self.matrix = self.stack.pop(0).matrix
        
class Axis3D(object):
    # construtor da Classe
    def __init__(self, x, y, z):
        self.v = np.array([x, y, z])
    
    def norm(self):
        return (self.v**2).sum()**0.5
        
    def getUnityVector(self):
        temp = self.v / self.norm()
        
        return Axis3D(temp[0], temp[1], temp[2])
    
    @staticmethod
    def crossProduct(a, b):
        res = np.cross(a.v, b.v)
        return Axis3D(res[0], res[1], res[2])
        
    
    @property
    def v(self):
        return self._v
    
    @v.setter
    def v(self, v):
        self._v = v
        
    def transform(self, T):
        temp = np.array([self.v[0], self.v[1], self.v[2], 1])
        print(temp)
        print(T.matrix)
        res = T.matrix.dot(temp.T)
        print(res)
        
        return Axis3D(res[0], res[1], res[2])

    def getPoint3D(self):
        return Point3D(self.v[0], self.v[1], self.v[2], 1, pygame.Color(0, 0, 0, 255), 1)
                      
class RotationInArbitraryAxis(CombinedTransformation):
    # construtor da classe
    def __init__(self, pointA, pointB, theta):
        super().__init__()
        super().add(Translation3D(-pointA.x, -pointA.y, -pointA.z))
        V = Axis3D(pointB.x - pointA.x, pointB.y - pointA.y, pointB.z - pointA.z)
        u = V.getUnityVector()
        d = np.sqrt(u.v[1]**2 + u.v[2]**2)
        alpha = np.arccos(u.v[2] / d)
        super().add(Rotation3D(alpha, 'x'))
        beta = np.arccos(d)
        super().add(Rotation3D(beta, 'y'))
        super().add(Rotation3D(theta, 'z'))
        super().add(Rotation3D(-beta, 'y'))
        super().add(Rotation3D(-alpha, 'x'))
        super().add(Translation3D(pointA.x, pointA.y, pointA.z))
        super().combine()


        
class ScaleWithFixedPoint(CombinedTransformation):
    # construtor da classe
    def __init__(self, point, Sx, Sy, Sz):
        super().__init__()
        super().add(Translation3D(-point.x, -point.y, -point.z))
        super().add(Scale3D(Sx, Sy, Sz))
        super().add(Translation3D(point.x, point.y, point.z))
        super().combine()
        
class WindowToViewportTransformation(CombinedTransformation):
    # construtor da classe
    def __init__(self, xd1, yd1, xd2, yd2, xv1, yv1, xv2, yv2):
        super().__init__()
        super().add(Translation3D(-xd1, -yd1, 0))
        x_view_factor = (xv2 - xv1) / (xd2 - xd1)
        y_view_factor = (yv2 - yv1) / (yd2 - yd1)
        super().add(Scale3D(x_view_factor, y_view_factor, 0))
        super().add(Translation3D(xv1, yv1, 0))
        super().combine()

class Picture2(object):
    # Construtor da Classe
    def __init__(self):
        self.primitivas = [] # Define uma lista de primitivas para representar um desenho
        self.count = 0
    
    def draw(self, screen):
        # Telhado
        p1 = Point3D(0, 350, 0, 1, pygame.Color(255, 0, 0, 255), 1)
        p2 = Point3D(700, 350, 0, 1, pygame.Color(255, 0, 0, 255), 2)
        p3 = Point3D(350, 0, 0, 1, pygame.Color(255, 0, 0, 255), 3)
        p4 = Point3D(350, 700, 0, 1, pygame.Color(255, 0, 0, 255), 4)
        l1 = Line(p1, p2, pygame.Color(255, 0, 0, 255), 1)
        l2 = Line(p3, p4, pygame.Color(255, 0, 0, 255), 2)
        pol = Polygon(True, pygame.Color(255, 0, 0, 255), True, pygame.Color(255, 255, 0, 255), 3)
        pol.addVertex(Point3D(450, 400, 0, 1, pygame.Color(0, 0, 0, 255), 1))
        pol.addVertex(Point3D(650, 450, 0, 1, pygame.Color(0, 0, 0, 255), 2))
        pol.addVertex(Point3D(575, 650, 0, 1, pygame.Color(0, 0, 0, 255), 3))
        
        T = RotationInArbitraryAxis(Point3D(350, 350, 0, 1, pygame.Color(255, 0, 0, 255), 1), Point3D(350, 350, 1, 1, pygame.Color(255, 0, 0, 255), 2), math.pi * self.count / 180)
        
        self.count = self.count + 10
        if self.count >= 360:
            self.count = 0
            
        pol2 = pol.transform(T)
        
        # Insere primitivas na lista
        self.primitivas.append(l1)
        self.primitivas.append(l2)
        self.primitivas.append(pol)
        self.primitivas.append(pol2)
        
        # Desenha cada primitiva que está na lista
        for item in self.primitivas:
            item.draw(screen)
        
        self.primitivas.clear()

import pygame

class Window(object):
    def __init__(self, xMin, yMin, xMax, yMax):
        self.xMin = xMin
        self.xMax = xMax
        self.yMin = yMin
        self.yMax = yMax
        
    @property
    def xMin(self):
        return self._xMin
    
    @xMin.setter
    def xMin(self, xMin):
        self._xMin = xMin
        
    @property
    def xMax(self):
        return self._xMax
    
    @xMax.setter
    def xMax(self, xMax):
        self._xMax = xMax
        
    @property
    def yMin(self):
        return self._yMin
    
    @yMin.setter
    def yMin(self, yMin):
        self._yMin = yMin
        
    @property
    def yMax(self):
        return self._yMax
    
    @yMax.setter
    def yMax(self, yMax):
        self._yMax = yMax
        
class Pipeline2D(object):
    def __init__(self):
        # Inicialização do PyGame
        pygame.init()            
        
    def run(self, listOfPrimitives, device, w, vp):
        while True:  # laço principal
            # captura eventos
            for event in pygame.event.get(): 
                
                # Captura evento de clicar em botão para fechar
                if event.type == pygame.QUIT:
                    return pygame.quit()
            
            l = []
            # Faz o recorte das primitivas
            for p in listOfPrimitives:
                l.extend(p.clip(window))
            
            # Mapeamento para Coordenadas Normalizadas
            T = WindowToViewportTransformation(w.xMin, w.yMin, w.xMax, w.yMax,-1.0, -1.0, 1.0, 1.0)
            #print(T.matrix)
            
            # Mapeamento para coordenadas de dispositivo
            T2 = WindowToViewportTransformation(-1.0, -1.0, 1.0, 1.0, vp.xMin, vp.yMin, vp.xMax, vp.yMax)
            #print(T2.matrix)
            T3 = CombinedTransformation()
            T3.add(T)
            T3.add(T2)
            T3.combine()
            #print(T3.matrix)
            
            # Rasterização e Exibição no dispositivo de saída
            # preencha a tela com a cor de fundo
            device.fill()
            
            #l = listOfPrimitives
            
            for p in l:
                p = p.transform(T3)
                p.draw(device)
            
            device.flip()

            pygame.display.update()
            
            pygame.time.Clock().tick(30)

class Picture3(object):
    # Construtor da Classe
    def __init__(self):
        self.primitivas = [] # Define uma lista de primitivas para representar um desenho
        self.count = 0
    
        # Telhado
        p1 = Point3D(-50, 300, 0, 1, pygame.Color(255, 0, 0, 255), 1)
        p2 = Point3D(100, 500, 0, 1, pygame.Color(255, 0, 0, 255), 2)
        p3 = Point3D(150, 600, 0, 1, pygame.Color(255, 0, 0, 255), 3)
        p4 = Point3D(250, 200, 0, 1, pygame.Color(255, 0, 0, 255), 4)
        p5 = Point3D(280, -50, 0, 1, pygame.Color(255, 0, 0, 255), 5)
        p6 = Point3D(750, 250, 0, 1, pygame.Color(255, 0, 0, 255), 6)
        p7 = Point3D(130, 800, 0, 1, pygame.Color(255, 0, 0, 255), 7)
        p8 = Point3D(800, 800, 0, 1, pygame.Color(255, 0, 0, 255), 8)
        p9 = Point3D(-100, 650, 0, 1, pygame.Color(255, 0, 0, 255), 9)
        p10 = Point3D(10, 800, 0, 1, pygame.Color(255, 0, 0, 255), 10)
        p11 = Point3D(0, 0, 0, 1, pygame.Color(255, 0, 0, 255), 11)
        p12 = Point3D(700, 0, 0, 1, pygame.Color(255, 0, 0, 255), 12)
        p13 = Point3D(700, 700, 0, 1, pygame.Color(255, 0, 0, 255), 13)
        p14 = Point3D(0, 700, 0, 1, pygame.Color(255, 0, 0, 255), 14)
        
        l1 = Line(p1, p2, pygame.Color(255, 0, 0, 255), 1)
        l2 = Line(p3, p4, pygame.Color(255, 0, 0, 255), 2)
        l3 = Line(p5, p6, pygame.Color(255, 0, 0, 255), 3)
        l4 = Line(p7, p8, pygame.Color(255, 0, 0, 255), 4)
        l5 = Line(p9, p10, pygame.Color(255, 0, 0, 255), 5)
        l6 = Line(p11, p12, pygame.Color(255, 0, 255, 255), 6)
        l7 = Line(p12, p13, pygame.Color(255, 0, 255, 255), 7)
        l8 = Line(p13, p14, pygame.Color(255, 0, 255, 255), 8)
        l9 = Line(p14, p11, pygame.Color(255, 0, 255, 255), 9)
        
        pol = Polygon(True, pygame.Color(255, 0, 0, 255), True, pygame.Color(255, 255, 0, 255), 3)
        pol.addVertex(Point3D(-50, 400, 0, 1, pygame.Color(0, 0, 0, 255), 1))
        pol.addVertex(Point3D(30, 20, 0, 1, pygame.Color(0, 0, 0, 255), 2))
        pol.addVertex(Point3D(450, -60, 0, 1, pygame.Color(0, 0, 0, 255), 3))
        pol.addVertex(Point3D(350, 730, 0, 1, pygame.Color(0, 0, 0, 255), 4))
        
        # Insere primitivas na lista
        self.primitivas.append(l1)
        self.primitivas.append(l2)
        self.primitivas.append(l3)
        self.primitivas.append(l4)
        self.primitivas.append(l5)
        self.primitivas.append(l6)
        self.primitivas.append(l7)
        self.primitivas.append(l8)
        self.primitivas.append(l9)
        self.primitivas.append(pol)
        
        # Desenha cada primitiva que está na lista
        #for item in self.primitivas:
        #    item.draw(screen)
        
        #self.primitivas.clear()
        
    def getListOfPrimitives(self):
        return self.primitivas

class Device(object):

    # Construtor da classe
    def __init__(self, title, bgColor, width, height):
        self.title = title       # título da janela
        self.bgColor = bgColor   # cor de fundo
        self.width = width       # largura da janela
        self.height = height     # altura da janela
        self.screen = pygame.display.set_mode(self.size()) # define o tamanho da tela
        pygame.display.set_caption(self.title)             # define o título da janela
        self.clock = pygame.time.Clock()
        
    def color(self):
        return self.bgColor
    
    def draw(self, primitive):
        primitive.draw(self)
        
    def fill(self):
        self.screen.fill(self.bgColor)
        
    def flip(self):
        self.screen.blit(pygame.transform.flip(self.screen, False, True), self.size())
        
    # retorna um vetor com o tamanho da tela
    def size(self):
        return (self.width, self.height)
    
    # modifica um pixel na tela com a cor desejada
    def setPixel(self, x, y, color):
        self.screen.set_at((x, y), color)
    
    # filtro da média para o antialising
    def meanFilter(self):
        # Captura a matrix da tela
        #frameBuffer2 = pygame.PixelArray(self.screen)
        
        from copy import copy
        frameBuffer = pygame.surfarray.array3d(self.screen)
        #print(frameBuffer)
        
        import numpy as np
        mask = np.ones((3, 3)) * 1/9 
        
        #print(mask)
               
        for i in range(1, self.width - 1):
            for j in range(1, self.height - 1):               
                temp = np.zeros((3))
                
                for k in range(-1,2):
                    for l in range(-1,2):
                        for b in range(3):
                            temp[b] = temp[b] + frameBuffer[i + k][j + l][b] * mask[k + 1][l + 1]
                        
                #print(pygame.Color(int(temp[0]), int(temp[1]), int(temp[2]), 255))
                        
                self.setPixel(i, j, pygame.Color(int(temp[0]), int(temp[1]), int(temp[2]), 255));
                #frameBuffer2[i][j] = pygame.Color(int(temp[0]), int(temp[1]), int(temp[2]), 255);

class PerspectiveTransformation(GeometricTransformation):
    def __init__(self, w):
        a = 2 * w.near / (w.xMax - w.xMin)
        b = (w.xMax + w.xMin) / (w.xMax - w.xMin)
        c = 2 * w.near / (w.yMax - w.yMin)
        d = (w.yMax + w.yMin) / (w.yMax - w.yMin)
        e = -(w.far + w.near) / (w.far - w.near)
        f = - 2 * w.far * w.near / (w.far - w.near)
        super().__init__(np.array([[a, 0, b, 0], [0, c, d, 0], [0, 0, e, f], [0, 0, -1, 0]]))
        
class OrthographicTransformation(GeometricTransformation):
    def __init__(self, w):
        a = 2/(w.xMax - w.xMin)
        b = -(w.xMax + w.xMin) / (w.xMax - w.xMin)
        c = 2/(w.yMax - w.yMin)
        d = -(w.yMax + w.yMin) / (w.yMax - w.yMin)
        e = -2 / (w.far - w.near)
        f = -(w.far + w.near) / (w.far - w.near)
        super().__init__(np.array([[a, 0, 0, b], [0, c, 0, d], [0, 0, e, f], [0, 0, 0, 1]]))

        
class SymmetricPerspectiveTransformation(GeometricTransformation):
    def __init__(self, fovy, aspect, w):
        fovy = fovy / 2
        f = 1.0 / math.tan(fovy * math.pi / 180.0)
        
        a = f / aspect
        b = f
        c = -(w.far + w.near) / (w.far - w.near)
        d = - 2 * w.far * w.near / (w.far - w.near)
        super().__init__(np.array([[a, 0, 0, 0], [0, b, 0, 0], [0, 0, c, d], [0, 0, -1, 0]]))
        
class Frustum(Window):
    def __init__(self, xMin, yMin, xMax, yMax, near, far):
        super().__init__(xMin, yMin, xMax, yMax)
        self.near = near
        self.far = far
    
    @property
    def near(self):
        return self._near
    
    @near.setter
    def near(self, near):
        self._near = near
        
    @property
    def far(self):
        return self._far
    
    @far.setter
    def far(self, far):
        self._far = far
          

class RotationInArbitraryPoint(CombinedTransformation):
    # construtor da classe
    def __init__(self, point, d, theta):
        super().__init__()
        super().add(Translation3D(-point.x, -point.y, -point.z))        
        super().add(Rotation3D(theta, d))
        super().add(Translation3D(point.x, point.y, point.z))
        super().combine()

class Camera(object):
    def __init__(self, eye, up, at):
        self.eye = eye # x0, y0, z0
        
        # sistema de coordenadas da câmera
        self.up = up
        self.at = at
        
        self.generateVisualizationTransform()
    
    def generateVisualizationTransform(self):
        self.T = CombinedTransformation()
        
        temp = Axis3D(self.eye.v[0] - self.at.v[0], self.eye.v[1] - self.at.v[1], self.eye.v[2] - self.at.v[2])
        #print(temp.v)
        #print(temp.norm())
        # Gera vetores unitários
        self.n = temp.getUnityVector()
        
        #print('Novo')
        #print(self.n.v)
        print(self.at.getUnityVector().v)
        
        temp = Axis3D.crossProduct(self.up, self.n)
        #print(temp.v)
        #normTemp = self.yv.norm()
        #print(normTemp)
        self.u = temp.getUnityVector()
        #print(self.u.v)
        #print(self.u.norm())
        
        self.v = Axis3D.crossProduct(self.n, self.u)
        #print(self.v.norm())
        
        #print(self.v.v)
        # Translada a câmera para origem
        self.T.add(Translation3D(-self.eye.v[0], -self.eye.v[1], -self.eye.v[2]))
        
        # Rotaciona o sistema de coordenadas cartesianas para o sistema da câmera
        #ORIGINAL
        self.T.add(GeometricTransformation(np.array([[self.u.v[0], self.u.v[1], self.u.v[2], 0], [self.v.v[0], self.v.v[1], self.v.v[2], 0], [self.n.v[0], self.n.v[1], self.n.v[2], 0], [0, 0, 0, 1]])))
        
        # Gera a matriz combinada
        self.T.combine()
        
    def getVisualizationTransform(self):
        return self.T
    
    def move_up(self):
        self.eye.v[1] += 10
        self.generateVisualizationTransform()
        
    def move_down(self):
        self.eye.v[1] -= 10
        self.generateVisualizationTransform()
        
    def move_right(self):
        self.eye.v[0] += 10
        self.generateVisualizationTransform()
    
    def move_left(self):
        self.eye.v[0] -= 10
        self.generateVisualizationTransform()
        
    def move_front(self):
        self.eye.v[2] += 10
        self.generateVisualizationTransform()
    
    def move_back(self):
        self.eye.v[2] -= 10
        self.generateVisualizationTransform()
        
    def rotate_right(self, angle):
        T = RotationInArbitraryPoint(self.at.getPoint3D(), 'y', angle)
        self.eye = self.eye.transform(T)
        self.generateVisualizationTransform()
        

class Pipeline3D(object):
    def __init__(self):
        # Inicialização do PyGame
        pygame.init()            
        self.count = 0
        
    def run(self, listOfPrimitives, device, w, vp, camera, proj):
        while True:  # laço principal
            # captura eventos
            for event in pygame.event.get(): 
                
                # Captura evento de clicar em botão para fechar
                if event.type == pygame.QUIT:
                    return pygame.quit()
                elif event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_UP:
                        print("UP")
                        camera.move_up()
                    elif event.key == pygame.K_DOWN:
                        print("DOWN")
                        camera.move_down()
                    elif event.key == pygame.K_RIGHT:
                        print("RIGHT")
                        camera.move_right()
                    elif event.key == pygame.K_LEFT:
                        print("LEFT")
                        camera.move_left()
                    elif event.key == pygame.K_PAGEUP:
                        print("FRONT")
                        camera.move_front()
                    elif event.key == pygame.K_PAGEDOWN:
                        print("BACK")
                        camera.move_back()
                    elif event.key == pygame.K_r:
                        print("RIGHT ROTATION")
                        self.count = self.count + 10
                        if self.count >= 360:
                            self.count = 0
                        camera.rotate_right(math.pi * self.count / 180)
                    
            TC = CombinedTransformation()
            
            # Transformação de Visualização
            T = camera.getVisualizationTransform()
            
            print(T.matrix)
            
            # Transformação de Projeção
            print('projection')
            print(proj.matrix)
            
            TC.add(T)
            TC.add(proj)
            TC.combine()
            
            print('TC')
            print(TC.matrix)
            
            l = []
            for p in listOfPrimitives:
                print('P')
                #print(p.getPoint())
                p2 = p.transform(TC)
                #print(p2.getPoint())
                p2.genRealPoints()
                #print(p2.getPoint())
                l.append(p2)
            
            # Mapeamento para coordenadas de dispositivo
            T = WindowToViewportTransformation(-1.0, -1.0, 1.0, 1.0, vp.xMin, vp.yMin, vp.xMax, vp.yMax)
            
            l2 = []
            # Faz o recorte das primitivas
            window = Window(vp.xMin, vp.yMin, vp.xMax, vp.yMax)
            for p in l:
                p2 = p.transform(T)
                l2.extend(p2.clip(window))
            
            
            
            # Rasterização e Exibição no dispositivo de saída
            # preencha a tela com a cor de fundo
            device.fill()
            
            #l = listOfPrimitives
            
            for p in l2:
                p.draw(device)
            
            # aplica o antialiasing
            # self.meanFilter()
            
            #device.flip()
            
            # atualiza a tela 
            pygame.display.update()
            
            pygame.time.Clock().tick(30)

class BezierCurve(Primitive):
    # construtor da classe
    def __init__(self, listOfControlPoints, step, color, id):
        super().__init__(id, 'bezier_curve')
        self.listOfControlPoints = listOfControlPoints # pontos de controle
        self.degree = 3       # assumindo cúbico
        self.color = color # cor do objeto
        self.M_B = GeometricTransformation(np.array([[-1, 3, -3, 1], [3, -6, 3, 0], [-3, 3, 0, 0], [1, 0, 0, 0]]))
        self.step = step
        
    
    # renderiza a curva desejada na tela
    def draw(self, screen):
        for u in np.linspace(0, 1, self.step):
            T = GeometricTransformation(np.array([[u**3, u**2, u, 1]]))
            comb = CombinedTransformation()
            comb.add(self.M_B)
            comb.add(T)
            comb.combine()
            print('Comb')
            print(comb.matrix)
            
            for i in range(0, len(self.listOfControlPoints)-1, 3):
                print(i)
                temp = np.array([[self.listOfControlPoints[i].x, self.listOfControlPoints[i].y, self.listOfControlPoints[i].z, self.listOfControlPoints[i].w], \
                                [self.listOfControlPoints[i + 1].x, self.listOfControlPoints[i + 1].y, self.listOfControlPoints[i + 1].z, self.listOfControlPoints[i + 1].w], \
                                [self.listOfControlPoints[i + 2].x, self.listOfControlPoints[i + 2].y, self.listOfControlPoints[i + 2].z, self.listOfControlPoints[i + 2].w], \
                                [self.listOfControlPoints[i + 3].x, self.listOfControlPoints[i + 3].y, self.listOfControlPoints[i + 3].z, self.listOfControlPoints[i + 3].w]])
                print(temp)
                
                res = comb.matrix.dot(temp)
                res = res.T
                print(res)
        
                screen.setPixel(int(res[0][0]), int(res[1][0]), self.color)           
            
    def transform(self, T):
        l = []
        
        for item in self.listOfControlPoints:
            point = item.getArray().T
            #print(T.matrix)
            res = T.matrix.dot(point)
            #print(res)
            l.append(Point3D(res[0], res[1], res[2], res[3], item.color, item.id))
        
        c = BezierCurve(l, self.step, self.color, self.id)
        return c
    
    def genRealPoints(self):
        for item in self.listOfControlPoints:
            item.divideByHomogeneousCoord()
        
    def clip(self, window): # TODO
        return [self]

class Picture5(object):
    # Construtor da Classe
    def __init__(self):
        self.primitivas = [] # Define uma lista de primitivas para representar um desenho
        self.count = 0
        
        
        # Pontos de Controle
        p1 = Point3D(0, 0, 20, 1.0, pygame.Color(255, 0, 0, 255), 2)
        p2 = Point3D(300, 300, 20, 1.0, pygame.Color(255, 0, 0, 255), 2)
        p3 = Point3D(700, 100, 20, 1.0, pygame.Color(255, 0, 0, 255), 2)
        p4 = Point3D(500, -100, 20, 1.0, pygame.Color(255, 0, 0, 255), 2)
        p5 = Point3D(400, 50, 20, 1.0, pygame.Color(255, 0, 0, 255), 2)
        p6 = Point3D(100, 100, 20, 1.0, pygame.Color(255, 0, 0, 255), 2)
        p7 = Point3D(0, 0, 20, 1.0, pygame.Color(255, 0, 0, 255), 2)
    
        l = []
        l.append(p1)
        l.append(p2)
        l.append(p3)
        l.append(p4)
        l.append(p5)
        l.append(p6)
        l.append(p7)
        
        curve = BezierCurve(l, 100, pygame.Color(0, 0, 0, 255), 3)
                
        # Insere primitivas na lista
        self.primitivas.append(curve)
        
    def getListOfPrimitives(self):
        return self.primitivas
    
screen = Device("Tela", pygame.Color(255, 255, 255, 255), 700, 700)
window = Frustum(-700, -700, 700, 700, 10.0, 100.0)
viewport = Window(0, 0, 700, 700)

# Para linhas
eye = Axis3D(350,350,-50)
up = Axis3D(0.0, 1.0, 0.0)
at = Axis3D(350, 350, 500)

cam = Camera(eye, up, at)

proj = OrthographicTransformation(window)
# Uso 1
#proj = SymmetricPerspectiveTransformation(170, 1, window)
# Uso 2 
#proj = SymmetricPerspectiveTransformation(60, 1, window)
#proj = PerspectiveTransformation(window)

pic = Picture5()

pip = Pipeline3D()

pip.run(pic.getListOfPrimitives(), screen, window, viewport, cam, proj)