#10.009 The Digital World 1D Project
#17F04 Group 2

''' MAIN SCRIPT '''
'''
This script takes in the adj_list from the 
subscriber script and the camera data, and 
activates the LEDs according to the various 
parameters. The usage is then published to 
be displayed on the Kivy app.
'''
#==============================================================================

'''============================================================================
                         Part I: Initialise program
============================================================================'''
import paho.mqtt.client as mqttClient #import the client1 for Google cloud
from collections import deque
import numpy as np
import argparse
import imutils
from tlkh_camera import PiVideoStream
import cv2
from time import sleep
#import RPi.GPIO as GPIO

'''=================================================== Setting up the camera'''
cv2.setUseOptimized(True) #OpenCV

#construct the argument parse and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-v", "--video",
    help="path to the (optional) video file")
ap.add_argument("-b", "--buffer", type=int, default=64,
    help="max buffer size")
args = vars(ap.parse_args())

#define the lower and upper boundaries of the "green"
#ball in the HSV color space, then initialize the
#list of tracked points
greenLower = (29, 86, 6)
greenUpper = (64, 255, 255)
pts = deque(maxlen=args["buffer"])

camera = PiVideoStream().start()
sleep(1.0)

def get_obj_position_and_brightness():
    # grab the current frame
    frame = camera.read()

    # resize the frame, blur it, and convert it to the HSV
    # color space
    frame = imutils.resize(frame, width=600)
    # blurred = cv2.GaussianBlur(frame, (11, 11), 0)
    hsv = cv2.cvtColor(frame, cv2.COLOR_BGR2HSV)
    #print(np.shape(hsv))
    #print(np.sum(np.sum(hsv[:,:,2])))
    # construct a mask for the color "green", then perform
    # a series of dilations and erosions to remove any small
    # blobs left in the mask
    mask = cv2.inRange(hsv, greenLower, greenUpper)
    mask = cv2.erode(mask, None, iterations=2)
    mask = cv2.dilate(mask, None, iterations=2)

    # find contours in the mask and initialize the current
    # (x, y) center of the ball
    cnts = cv2.findContours(mask.copy(), cv2.RETR_EXTERNAL,
        cv2.CHAIN_APPROX_SIMPLE)[-2]
    center = None
    
    x,y = 0,0

    # only proceed if at least one contour was found
    if len(cnts) > 0:
        # find the largest contour in the mask, then use
        # it to compute the minimum enclosing circle and
        # centroid
        c = max(cnts, key=cv2.contourArea)
        ((x, y), radius) = cv2.minEnclosingCircle(c)
        M = cv2.moments(c)
        center = (int(M["m10"] / M["m00"]), int(M["m01"] / M["m00"]))

        # only proceed if the radius meets a minimum size
        if radius > 10:
            # draw the circle and centroid on the frame,
            # then update the list of tracked points
            cv2.circle(frame, (int(x), int(y)), int(radius),
                (0, 255, 255), 2)
            cv2.circle(frame, center, 5, (0, 0, 255), -1)
            
            general_brightness = np.sum(np.sum(hsv[:,:,2]))/345600.
            sleep(0.1)
        
        # update the points queue
        pts.appendleft(center)
        # loop over the set of tracked points
        for i in xrange(1, len(pts)):
            # if either of the tracked points are None, ignore
            # them
            if pts[i - 1] is None or pts[i] is None:
                continue

            # otherwise, compute the thickness of the line and
            # draw the connecting lines
            thickness = int(np.sqrt(args["buffer"] / float(i + 1)) * 2.5)
            cv2.line(frame, pts[i - 1], pts[i], (0, 0, 255), thickness)

        # show the frame to our screen
        cv2.imshow("Frame", frame)
    
    person = Target(int(x), int(y))
    print('Person is at {}; brightness is {}'.format(person, general_brightness))
    return person, general_brightness

'''============================================== Setting up the LED control'''
#set up coordinate system
class Target: #for person
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y
    def __str__(self):
        return '({},{})'.format(self.x, self.y)

class LED_coord: #for LEDs
    def __init__(self, x, y, pin, userpref = 0.8):
        self.x = x
        self.y = y
        self.pin = pin           #GPIO pin number
        self.userpref = userpref #user preference of LED's brightness (0-1)
        self.dist = 0.           #distance from person
        self.duty = 100          #duty cycle varying with brightness
    def __str__(self):
        return 'LED ({},{}), GPIO pin {}'.format(self.x, self.y, self.pin)

#set LED coordinates
p1 = LED_coord(120,112,19)         #   legend  # - LED
p2 = LED_coord(240,112,5)          #           x - camera (600 x 450 px)
p3 = LED_coord(360,112,13)         #   origin O-------------------+ 0
p4 = LED_coord(480,112,26)         #          |   |   |   |   |   |
p5 = LED_coord(120,225,12)         #          |--1#--2#--3#--4#---| 112
p6 = LED_coord(240,225,21)         #          |   |   |   |   |   |
p7 = LED_coord(360,225,20)         #          |--5#--6#-x7#--8#---| 225
p8 = LED_coord(480,225,16)         #          |   |   |   |   |   |
p9 = LED_coord(120,338,6)          #          |--9#-10#-11#-12#---| 338
p10 = LED_coord(240,338,7)         #          |   |   |   |   |   |
p11 = LED_coord(360,338,8)         #          +-------------------+ 450
p12 = LED_coord(480,338,9)         #          0  120 240 360 480 600
led_list = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12]

'''======================================================= Reading text file'''
#directory of text file receiving adj_list data
txt = 'txt'

def get_adj_list():
    f = open(txt, 'r')
    reading = f.readline()
    print('Received adj_list: {}'.format(reading))
    adj_list = reading.split(',')
    f.close()
    return adj_list

'''============================================================================
                           Part II: LED Operation
============================================================================'''

#decide brightness of each LED
def decide_brightness():
    #obtain person's position, brightness, adj_list
    person, brightness = get_obj_position_and_brightness()
    adj_list = get_adj_list()
    
    #alter base attribute of LEDs
    for i in range(len(led_list)):
        led_list[i].userpref = adj_list[i]
    #for debugging=============================================================
    check_userpref = [led.userpref for led in led_list]
    print('LED adjustments: {}'.format(check_userpref))
    #==========================================================================
    
    #determine duty cycle of LEDs
    check_usage = ''
    for led in led_list:
        #alter dist attribute
        led.dist = ( (person.x - led.x)**2 + (person.y - led.y)**2 )**0.5
        #alter duty attribute
        power1 = (100 - brightness) * led.userpref  #factoring user pref
        if led.dist < 165:    #--------------LED is at most diagonally adjacent
            mult = ( (165 - led.dist) / 165 )**0.25
            led.duty = power1 * mult
            check_usage += '{},'.format(led.duty / 100)
        else:    #---------------------------LED is too far away, do not use
            led.duty = 0
            check_usage += '0,'
        check_usage.rstrip(',')
    
    #for publishing and debugging==============================================
    print('LED usage: {}'.format(check_usage))
    dw1d.publish("duty_list", str(check_usage))
    #==========================================================================

#activate LEDs according to duty cycle assigned
def activate_led():
    GPIO.setmode(GPIO.BCM)
    print('Activating LEDs')
    for led in led_list:
        GPIO.setup(led.pin, GPIO.OUT)
        pwm_led = GPIO.PWM(led.pin, 50)
        pwm_led.start(100)
        pwm_led.ChangeDutyCycle(led.duty)
    GPIO.cleanup()

'''============================================================================
           Part III: Actual Operation & Publishing to Google Cloud
============================================================================'''
def on_connect(client, userdata, flags, rc):
    if rc == 0:
         print("Connected to broker")
         global Connected
         Connected = True
    else:
         print("Connection failed")

#setting up connection to Google Cloud
Connected = False
broker_address="35.197.131.13"
port = 8883
print("Creating new instance")
dw1d = mqttClient.Client("DW1D")
dw1d.username_pw_set("sammy","password")  #set usernames and passwords
dw1d.on_connect = on_connect
print("Connecting to broker")
dw1d.connect(broker_address, port=port)   #connect to broker

while True:
    decide_brightness()
    activate_led()

#cleanup
print("Stopping camera")
cv2.destroyAllWindows()
camera.stop()
print("Cleaning GPIO")
GPIO.cleanup()
print("=== END ===")