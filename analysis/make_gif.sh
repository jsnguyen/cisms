#!/bin/bash

echo "Making gif..."
ffmpeg -y -r 24 -i ./frames/%05d.png ./animation.gif
