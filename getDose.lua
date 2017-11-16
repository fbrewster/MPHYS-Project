--Check Heart Calcifications
--[[
Frank Brewster
Finds calcifications and then checks what percentage of these are within the heart delineation
]]--


--basefolder = [[D:\MPHYS_Data\]]
basefolder = [[D:\MPHYS\Data\]]
currentpatientpack = [[200509064.pack]]
clipDist = 0.5--acceptabel distance from the centre of a bubble to the CH boundary
maxVol = 5--largest volumes considered a calcification
bleedVol = 10--volume above which a fill is considred to have bled
maxVolTot = 0--to keep track of how big the bleed volumes are

function getRandScan()--Gets a random scan name from the list in ..\Data
  local file = io.open(basefolder..[[list.txt]])--open list
  local fNames = {}
  local i = 1

  if file then--if the file exists
    for line in file:lines() do--iterate through the lines
        fNames[i] = unpack(line:split(" "))--split by line
        i=i+1
    end
  end
  
  local nOfFiles = i-1--number of files to choose from
  math.randomseed(os.time())--initiate a random number generator
  local dummy = math.random(nOfFiles)--first number can be non-random
  local index = math.random(nOfFiles)--get a random index between 1 and the number of files
  local randScan = fNames[index]--get the file name at that index
  
  return randScan
end

function TScan:setup()
  self.adjust = mAdjust
  self.Transform = mTrans
end

function display(scans)
  for i=1,#scans do
    local j= i+2
    wm.Scan[j] = scans[i]:copy()
  end
  newlinkedpage:activate()
end


function lungConvexHull()--Makes a convex hull from the delineation of both lungs
  
  local lungs = wm.Delineation.Both_Lungs or wm.Delineation.Lungs or wm.Delineation.Both_Lung --get lung delineation
   
  local chIndex = Field:new()--create index for CH
  AVS:CONVEX_HULL( lungs.Dots, chIndex )--make CH index from lung delin
  chDots = lungs.Dots--dots are the same as lung
  
  local trans = mTrans--take scan transform
  AVS:DOTXFM( chDots, trans, chDots )--change chDots to scan transform
  local delChPts = Field:new()
  local delChIndex = Field:new()
  AVS:TT_TO_CNT( chDots, chIndex, wm.Scan[1].Data, delChPts, delChIndex )--makes delin from dots and index
  AVS:TRANSFORM_MATH( trans, nil, trans, true )--inverts trans
  AVS:DOTXFM( delChPts, trans, delChPts )--applied inverted trans to get back
  AVS:VRML_WRITE( delChPts, delChIndex, nil, [[c:\temp\lungsConvexHull_contoursScan.wrl]] )--writes 3d model of ch to file
  
  --clean and fix contours
  results = String:new()
  AVS:CNT_CLEAN_VERTICES( delChPts, delChIndex, delChPts, delChIndex, results )
  AVS:CNT_CLEAN( delChPts, delChIndex, delChPts, delChIndex, true )
  
  local newDel = lungs:copy()
  newDel.Dots = delChPts
  newDel.Index = delChIndex
  
  local cHull = Scan:new()
  cHull.Adjust = mAdjust
  cHull.Transform = mTrans
  
  cHull = wm.Scan[1]:burn(newDel, 255, false)--burn delineation onto scan 3
  cHull.Description = "Lung Convex Hull Delineation"
  
  return newDel, cHull
end


function mask(input, cHull, lShrink)--masks scan[input] with scan[cHull], shrunk by "shrink" and outputs to scan[output]
  local ch = cHull:copy()
  lShrink = (lShrink*100)+127--convert to threshold value
  --Exclude data outside of CH
  local output = Scan:new() 
  output:setup()
  --need to shrink scan
  AVS:FIELD_OPS(ch.Data, ch.Data, AVS.FIELD_OPS_SignedDist)--make singed distance image
  local temp = ch:copy()
  AVS:FIELD_THRESHOLD(temp.Data, ch.Data, lShrink)--delete everything shrink cm in
  AVS:FIELD_MASK(input.Data, ch.Data, output.Data)--delete outside of new shrunk mask
  ch.Description = [[Convex hull shrunk by <shrink>]]
  output.Description = "Masked and thresholded data"
  return output
end



function chDist(bubble, clipDist, chIn)--Checks if a bubble is too close to the convex hull in scan 3

  
  AVS:FIELD_OPS(chIn.Data, chIn.Data, AVS.FIELD_OPS_SignedDist)--turn ch into a distance measure
  local centre = Field:new()
  AVS:FIELD_CENTER_DOT(bubble, centre)-- put the centre of bubble into centre
  local distF = Field:new()
  AVS:FIELD_SAMPLE(centre, chIn.Data, distF)-- sample the field ch at all the points in centre and output to dist
  --print(distF:getvalue(0))
  local dist = distF:getvalue(0)--extract distance value
  local distCent = (dist.value - 127)/100--convert to cm from CH boundary
  local distAbs = math.abs(distCent)--take absolute value to include inside and out
  --print(distAbs)
  
  local clip = true
  
  if distAbs<clipDist then-- if the distance from the centre of the bubble to the ch is less than clipDist the return false
     clip = false
  end
  
  return clip
end

function bubbleCent (inSc,cHull) 
  
  local dummy = field:new()
  local labels = Scan:new()
  labels:setup()
  local temp = Field:new()
  print("start bubble")
  AVS:FIELD_TO_INT( inSc.Data, labels.Data)
  AVS:FIELD_LABEL( labels.Data, labels.Data, dummy, AVS.FIELD_LABEL_3D, 1 )
  
  local  cents = {}
  local nBubbles = labels.Data:max().value
  local volTemp = Double:new()
  if nBubbles>2000 then
    wm.Scan[4] = inSc:copy()
    print("Too many bubbles! " .. nBubbles .. " Pack: " .. currentpatientpack)
    AVS:FIELD_OPS( inSc.Data, inSc.Data, 1, AVS.FIELD_OPS_Smooth )
    AVS:FIELD_TO_INT( inSc.Data, labels.Data)
    AVS:FIELD_LABEL( labels.Data, labels.Data, dummy, AVS.FIELD_LABEL_3D, 1 )
    nBubbles = labels.Data:max().value
    if nBubbles>2000 then
      wm.Scan[4] = inSc:copy()
      error("Tried smoothing again but still too many bubbles. " .. nBubbles)
    else
      print("Smoothing again worked")
    end
  end
  print('# of bubbles: ' .. nBubbles)

  
  for i = 1,nBubbles do
    volTemp = labels:volume(i,i).value
    if volTemp < maxVol and volTemp>0.001 then
      AVS:FIELD_THRESHOLD( labels.Data, temp, i, i, 255, 0)
      --local clip = chDist(temp,clipDist,cHull)
      if true then
        local tempCent = field:new()
        AVS:FIELD_CENTER_DOT(temp, tempCent)
        table.insert( cents, {id=i, cent=tempCent})
      end
    end
    if i == (nBubbles-nBubbles%2)/2 then
      print("Halfway through bubble checking")
    end
    collectgarbage()
  end
  
  local output = {}
  local outTot = Scan:new()
  local tempBubScan = Scan:new()
  tempBubScan:setup()
  outTot:setup()
  
  for i = 1, #cents do
    AVS:FIELD_THRESHOLD(labels.Data, temp, cents[i].id, cents[i].id, 255, 0)
    outTot.Data:add(temp)
    tempBubScan.data=temp
    output[i] = tempBubScan:copy()
  end
  
  return cents, outTot, output
end

function hasBled(inScan)
  local inThr = Scan:new()
  local labels = Scan:new()
  inThr.Adjust = mAdjust
  inThr.Transform = mTrans
  labels.Adjust = mAdjust
  labels.Transform = mTrans
  local dummy = field:new()
  local volTemp,volTot = Double:new()
  volTot = 0
  
  AVS:FIELD_THRESHOLD(inScan.Data, inThr.Data, 5000, 5000, 255, 0)
  
  AVS:FIELD_TO_INT(inThr.Data, labels.Data)
  AVS:FIELD_LABEL( labels.Data, labels.Data, dummy, AVS.FIELD_LABEL_3D, 1 )
  
  local nBubbles = labels.Data:max().value
  
  for i=1,nBubbles do
    volTemp = labels:volume(i,i).value
    volTot = volTot+volTemp
  end

  
  local bleed = false
  if volTot>=bleedVol then
    bleed = true
  elseif volTot>maxVolTot then
    maxVolTot=volTot
  end
  
  return bleed
end

--currentpatientpack = getRandScan()
loadpack(basefolder ..[[PackWithHearts\]] .. currentpatientpack)
--loadpack(basefolder .. [[Pack\]] .. currentpatientpack)
--wm.Scan[1].Description = currentpatientpack
local scans = {}

--set master adjust and transform.
mAdjust = wm.Scan[1].Adjust
mTrans = wm.Scan[1].Transform

local orign = Scan:new()
orign = wm.Scan[1]:copy()

--Make convex hull
local lungCH = Delineation:new()
local cHull = Scan:new()
cHull:setup()
lungCH,cHull = lungConvexHull()

local masked = Scan:new()
masked:setup()
masked = mask(orign,cHull,1)
local hist = masked:histogram(masked,650,3000,3000,256)
hist.cumulative = true
local lowThresh = hist:percentile(99)
--[[local floodThresh = hist:percentile(97)
if floodThresh.value<1135 then 
  print("Flood threshold's fucked mate")
  floodThresh = 1135
end]]
local highThresh = hist:max()
print("got hist")
--threshold into scan 4
local threshed = Scan:new()
threshed:setup()
AVS:FIELD_THRESHOLD( orign.Data, threshed.Data, lowThresh.value, highThresh.value )
scans[2] = threshed

--Smooth and put in scan 5
local smoothed = Scan:new()
smoothed:setup()
AVS:FIELD_OPS( threshed.Data, smoothed.Data, 2, AVS.FIELD_OPS_Smooth )--Smoothing filter
--local tempScan = Field:new()
local tempScan = smoothed.Data:copy()
AVS:FIELD_THRESHOLD( tempScan, smoothed.Data, 1, 255)--Flatten
AVS:FIELD_OPS(smoothed.Data, smoothed.Data, 5, AVS.FIELD_OPS_Closing)--Closing filter to minimise # of volumes
print("smoothed")
--spine burn
spine = Scan:new()
spine:setup()
local delin = wm.Delineation.SC or wm.Delineation.cord or wm.Delineation.SCanal or wm.Delineation.SCANAL or wm.Delineation.SCord or wm.Delineation.sc
spine = orign:burn(delin, 255, true)
AVS:FIELD_OPS(spine.Data, spine.Data, AVS.FIELD_OPS_SignedDist)
local spineBig = Scan:new()
AVS:FIELD_THRESHOLD(spine.Data, spineBig.Data,1)
smoothed.Data:add(spineBig.Data)
scans[3] = smoothed

--Mask with CH and put in scan 6. Shrink by 1cm
local maskSmooth = Scan:new()
maskSmooth:setup()
maskSmooth = mask(smoothed,cHull,1)
scans[4] = maskSmooth
print("masked")


local cents = field:new()
local bubTot = Scan:new()
bubTot:setup()
bubs = {}
local shcHull = cHull:copy()
AVS:FIELD_OPS(shcHull.Data, shcHull.Data, AVS.FIELD_OPS_SignedDist)--make singed distance image
tempScan = shcHull.Data:copy()
AVS:FIELD_THRESHOLD(tempScan, shcHull.Data, 227)--delete everything shrink cm in
--shcHull.Data = tempScan
scans[1] = shcHull
cents,bubTot,bubs = bubbleCent(maskSmooth,shcHull:copy())--find bubbles and centres
scans[5] = bubTot
print("Found centres, starting flooding")
local halfway = (#cents-#cents%2)/2
local flood = orign:copy()
local allMask = Scan:new()
allMask:setup()
local shitThresh = {}
--mask(8,3,8)

for i=1,#cents do--floodfill from centre of each bubble
  local currentBub = bubs[i]:copy()
  local bubMask = Scan:new()
  bubMask:setup()
  AVS:FIELD_OPS(currentBub.Data, currentBub.Data, 4, AVS.FIELD_OPS_Smooth3)
  AVS:FIELD_MASK(orign.Data, currentBub.Data, bubMask.Data)
  local bubHist = bubMask:histogram(bubMask,650,3000,3000,256)
  bubHist.cumulative = true
  local floodThresh = bubHist:percentile(70)
  --if floodThresh.value<1100 then 
    --floodThresh.value = 1135 
  --end
  allMask:add(bubMask)
  
  local safe = flood:copy()
  AVS:FLOODFILL(cents[i].cent, flood.Data, floodThresh.value, 3, 5000)
  local bleed = hasBled(flood)--check for bleed
  collectgarbage()
  if bleed then--if it has bled, don't include
    flood = safe:copy()
    shitThresh[i] = floodThresh.value 
    AVS:FLOODFILL(cents[i].cent, flood.Data, 1135, 3, 5000)
    local bleed2 = hasBled(flood)
    if bleed2 then
      flood = safe:copy()
    end
  end
  if i == halfway then
    print("Halway through flooding")
  end
end
scans[6] = flood
print("Finished flooding")

scans[1] = allMask
display(scans)

--[[
wm.Scan[1].Description = "Original CT Scan"
wm.Scan[2].Description = "Dose"
wm.Scan[3].Description = "Mask"
wm.Scan[4].Description = "Thresholded"
wm.Scan[5].Description = "Smoothed"
wm.Scan[6].Description = "Masked"
wm.Scan[7].Description = "Calcifications"
wm.Scan[8].Description = "Original with calcifications flooded to 5000"]]

local dummy = field:new()
--[[local bubFlood = Scan:new()
bubFlood:setup()
totBubFlood = Scan:new()
totBubFlood:setup()
AVS:FIELD_THRESHOLD( flood.Data, totBubFlood.Data, 4999,5001)
AVS:FIELD_OPS(totBubFlood.Data, totBubFlood.Data, 2, AVS.FIELD_OPS_Smooth)
tempScan = totBubFlood.Data:copy()
AVS:FIELD_THRESHOLD(tempScan, totBubFlood.Data, 1, 255)
AVS:FIELD_TO_INT( totBubFlood.Data, bubFlood.Data )
AVS:FIELD_LABEL( bubFlood.Data, bubFlood.Data, dummy, AVS.FIELD_LABEL_3D, 1)
local doseScan = wm.Scan[2]:copy()
local temp = Scan:new()
temp:setup()
local medianDose = {}
local maxDose = doseScan.Data:max().value
maxDose = (maxDose - maxDose%1) + 1
for i=1,#cents do
  AVS:FIELD_THRESHOLD( bubFlood.Data, temp.Data, i, i)
  local bubHist = temp:histogram(doseScan, 254, 256, 100, 50)
  local median = bubHist:percentile(50).value
  if median>0 then
    table.insert(medianDose, median)
  end
  if median == 1 then
    wm.Marker[1]:fit(temp)
    wm.Marker[1]:gotomarker()
    print("wierd median")
  end
  collectgarbage()
end]]


print("Stop")

--[[inHeartTot = inHeartTot + currentInHeart
bubblesTot = bubblesTot + #cents]]