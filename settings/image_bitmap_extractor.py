import numpy as np
import matplotlib.pyplot as plt
import cv2
from PIL import Image, ImageChops



def main():
    #Cargar Imagen
    image = cv2.imread('test1.jpg')
    plt.imshow(image)
    plt.show()

    #Para cortar la imagen
    '''
    cut_x0 = 300
    cut_xf = 2850
    cut_y0 = 700
    cut_yf = 3280
    image = image[cut_y0:cut_yf, cut_x0:cut_xf]
    plt.imshow(image)
    plt.show()
    '''

    #Extraer contorno
    image = cv2.cvtColor(image,cv2.COLOR_BGR2GRAY)
    image = cv2.GaussianBlur(image,(3,3), 2*cv2.BORDER_DEFAULT)
    image = cv2.convertScaleAbs(image, alpha=2.6, beta=50) 
    image = cv2.equalizeHist(image)

    flag, image = cv2.threshold(image, 40, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
                                       #Mover ese 40 según sea necesario

    #Quitar puntitos
    #image = cv2.GaussianBlur(image,(9,9), cv2.BORDER_DEFAULT)
    #flag, image = cv2.threshold(image, 254, 255, cv2.THRESH_BINARY)

    plt.imshow(image)
    plt.show()


    #Guardar imagen
    cv2.imwrite('threshold.bmp', image)

if __name__ == '__main__':
    main()