--Iterate through data

wm:clear()
basefolder = [[D:\MPHYS_Data\]]
nOfScans=30
inHeartTot=0
bubblesTot=0

local file = io.open(basefolder..[[list.txt]])--open list
local fNames = {}
local i = 1

if file then--if the file exists
  for line in file:lines() do--iterate through the lines
    fNames[i] = unpack(line:split(" "))--split by line
    i=i+1
  end
end
io.close(file)

file =io.open(basefolder..[[inHeart.csv]],"w")
file:write("id,# of cal,Frac of cal in heart \n")
io.close(file)

local nOfFiles = i-1

local indices = {}
math.randomseed(os.time())
local dummy = math.random(nOfFiles)
for j=1,nOfScans do
  indices[j]=math.random(nOfFiles)
end

local checkHeartCal = assert(loadfile([[C:\Users\Frank\MPHYS\code\MPHYS-Project\checkHeartCal.lua]]))
for j=1,nOfScans do
  local index = indices[j]
  currentpatientpack = fNames[index]
  print(j)
  checkHeartCal()
  wm:clear()
end

local heartProp = (inHeartTot/bubblesTot)*100
print(heartProp .. "% of calcifications found are in the heart")