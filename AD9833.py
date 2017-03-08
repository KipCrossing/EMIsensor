# main.py -- put your code here!


from pyb import Pin
from pyb import SPI
from pyb import ADC
import time


Freq =10000
(ClockFreq,x,x,x) = pyb.freq()


word = hex(round((Freq*2**28)/ClockFreq))


MSB = (int(word) & 0xFFFC000)>>14

LSB = int(word) & 0x3FFF

LSB |= 0x4000
MSB |= 0x4000


p_out = Pin('X5', Pin.OUT_PP)


spi = SPI(1, SPI.MASTER, baudrate=9600, polarity=1, phase=0,firstbit=SPI.MSB)


def bytes(integer):
    return divmod(integer, 0x100)


def Send(data):
    high, low = bytes(data)
    p_out.low()

    spi.send(high)
    spi.send(low)

    p_out.high()



Send(0x2100)
Send(LSB)
Send(MSB)
Send(0xC000)
Send(0x2000)    #0x2020 - square, 0x2000 - sin, 0x2002 - triangle






adc = ADC(Pin('X19'))
while True:
    print(adc.read()) # read value, 0-4095
