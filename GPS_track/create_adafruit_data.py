from Adafruit_IO import Client, Data
from pylab import *
import random
import time

aio2 = Client(username='jozsamate',key='aio_tIbu493l7whUSGyXPAsFNoJSD8kW')
# aio = Client(username='Ballon1',key='aio_ETIU18iWTGajN4wD2fawgjLbgwLi')
# Data = aio.data('ballon', max_results=10)
# data = Data[0].value

with open('sensor_data.txt', 'r') as file:
    data = file.read()

latitude, longitude = 46.7675, 23.5895

for i in range(10000):
	latitude  += 0.0002#(random.randint(0, 1)*2-1)*0.0001 + 0.0001
	longitude += 0.0002#(random.randint(0, 1)*2-1)*0.0001
	print(latitude, longitude)
	data2 = data.split("Lati= ")[0]+"Lati= "+str(latitude)+", Longi= "+str(longitude)+", Altitude= "+data.split("Altitude= ")[1]
	Data2 = Data(value=data2)
	aio2.create_data("randomwalk2", Data2)
	time.sleep(3)