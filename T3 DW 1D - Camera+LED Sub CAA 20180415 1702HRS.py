#10.009 The Digital World 1D Project
#17F04 Group 2

''' SUBSCRIBER SCRIPT'''
'''
This script subscribes to the adj_list topic
on Google Cloud. The received input is written
to a text file, which is then read by the main
script to perform the LED operations.
'''
#directory of text file receiving adj_list data
txt = 'txt'

#import statements
import paho.mqtt.client as mqttClient #import the client1 for Google cloud
from time import sleep

#write received data in a text file
def write_to_txt(file, message):
    f = open(file, 'w')
    f.write('{}\n'.format(message))
    f.close()
    print('Message written to text file')

#call write_to_txt() upon receiving message from Google Cloud
def on_message(client, userdata, message):
    received_data = str(message.payload.decode("utf-8"))
    print("Message received:", received_data)
    write_to_txt(txt, received_data)

#setting up connection to Google Cloud
broker_address="35.197.131.13"
port = 8883
print("Creating new instance")
dw1d = mqttClient.Client("DW1D")
dw1d.username_pw_set("sammy","password")  #set usernames and passwords
dw1d.on_message = on_message              #attach functions to callback
print("Connecting to broker")
dw1d.connect(broker_address, port=port)   #connect to broker

#actual loop for receiving info
dw1d.loop_start()
dw1d.subscribe("adj_list")
sleep(100)
dw1d.loop_stop()