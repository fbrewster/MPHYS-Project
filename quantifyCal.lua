--Quantify Calcifications
--[[
Frank Brewster
Creates convex hull around both lungs delineation and then thresholds and masks.
Then does connected pixel analysis and discards large volumes.
Floods volumes to minimum calcification CT number.
]]--

local testpatientpack = [[201110978.pack]]--scan to use if testing a specific scan
basefolder = [[C:\Users\Frank\MPHYS\Data\]]--where the data is
clipDist = 0.3--acceptabel distance from the centre of a bubble to the CH boundary
maxVol = 2.7--largest volumes considered a calcification
bleedVol = 75--volume above which a fill is considred to have bled
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

function lungConvexHull()--Makes a convex hull from the delineation of both lungs
  
  local lungs = wm.Delineation.Both_Lungs or wm.Delineation.Lungs or wm.Delineation.Both_Lung --get lung delineation
   
  local chIndex = Field:new()--create index for CH
  AVS:CONVEX_HULL( lungs.Dots, chIndex )--make CH index from lung delin
  chDots = lungs.Dots--dots are the same as lung
  
  trans = wm.Scan[1].Transform--take scan transform
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
  
  wm.Scan[3] = wm.Scan[1]:burn(newDel, 255, false)--burn delineation onto scan 3
  wm.Scan[3].Description = "Lung Convex Hull Delineation"
  
  return newDel
end


function mask(input, cHull, output, shrink)--masks scan[input] with scan[cHull], shrunk by "shrink" and outputs to scan[output]
  if shrink then
    shrink = (shrink*100)+127--convert to threshold value
  else
    local shrink = 0
  end
  --Exclude data outside of CH
  wm.Scan[output].Adjust = wm.Scan[input].Adjust
  wm.Scan[output].Transform = wm.Scan[input].Transform
  --need to shrink scan
  AVS:FIELD_OPS(wm.Scan[cHull].Data, wm.Scan[cHull].Data, AVS.FIELD_OPS_SignedDist)--make singed distance image
  AVS:FIELD_THRESHOLD(wm.Scan[cHull].Data, wm.Scan[cHull].Data, shrink)--delete everything shrink cm in
  AVS:FIELD_MASK(wm.Scan[input].Data, wm.Scan[cHull].Data, wm.Scan[output].Data)--delete outside of new shrunk mask
  wm.Scan[cHull].Description = [[Convex hull shrunk by <shrink>]]
  wm.Scan[output].Description = "Masked and thresholded data"
end

--Messy obsolete function for excluding large volumes
function excludeVols(inSc, outSc)
  wm.Scan[outSc].Adjust = wm.Scan[inSc].Adjust
  wm.Scan[outSc].Transform = wm.Scan[inSc].Transform
  local dummy = field:new()
  local labels = Scan:new()
  labels.Adjust = wm.Scan[inSc].Adjust
  labels.Transform = wm.Scan[inSc].Transform
  local temp = Field:new()
  
  AVS:FIELD_TO_INT( wm.Scan[inSc].Data, labels.Data)
  AVS:FIELD_LABEL( labels.Data, labels.Data, dummy, AVS.FIELD_LABEL_3D, 1 )
  
  local vols = {}
  local nBubbles = labels.Data:max().value
  if nBubbles>2000 then
    error("Too many bubbles! " .. nBubbles)
  end
  print('# of bubbles: ' .. nBubbles)
  local volTemp
  
  local file = io.open(basefolder..[[bubbles.csv]],"w")
  if file then
    file:write("id,vol,x,y,z \n")
  else
    error([[Error opening bubbles file, is it already open?]])
  end
  
  for i = 1,nBubbles do
    volTemp = labels:volume(i,i).value
    AVS:FIELD_THRESHOLD( labels.Data, temp, i, i, 255, 0)
    local x,y,z = Double:new()
    x,y,z = temp:center()
    file:write(i,',',volTemp,',',x,',',y,',',z,'\n' )
    if volTemp < maxVol and volTemp>0.01 then
      local clip = chDist(temp,clipDist)
      --print(clip)
      if clip then
        table.insert( vols, { id=i, vol=volTemp } )
      end
    end
  end
  io.close(file)
  
  
  local function compareVols(a,b)
    return a.vol > b.vol
  end
  table.sort(vols, compareVols)
  
  

  for j=1, #vols do
    AVS:FIELD_THRESHOLD( labels.Data, temp, vols[j].id, vols[j].id, 255, 0)
    --if j==1 then
      --local centre = temp:center()
      --print("Largest bubble located at " .. centre)
    --[[local x,y,z = Double:new()
    x,y,z = temp:center()
    wm.Scan[outSc].Data:add(temp)
    file:write(vols[j].id,',',vols[j].vol,',',x,',',y,',',z,'\n' )]]
    wm.Scan[outSc].Data:add(temp)
  end
  --io.close(file)
  wm.Scan[6].Description = "Volumes over 0.01cm^3 and under maxVol cm^3"
  print("# of final vols: " .. #vols)
  if #vols~=0 then
    print("Largest vol: " .. vols[1].vol)
  else
    error("No volumes detected")
  end
end

function chDist(bubble, clipDist)--Checks if a bubble is too close to the convex hull in scan 3
  local ch = wm.Scan[3]:copy()
  local chIn = wm.Scan[3]:copy()
  local chOut = wm.Scan[3]:copy()

  
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

function bubbleCent (inSc)--Finds centres of bubbles
  
  local dummy = field:new()
  local labels = Scan:new()
  labels.Adjust = mAdjust
  labels.Transform = mTrans
  local temp = Field:new()
  
  AVS:FIELD_TO_INT( wm.Scan[inSc].Data, labels.Data)--Changes pixel values to integers
  AVS:FIELD_LABEL( labels.Data, labels.Data, dummy, AVS.FIELD_LABEL_3D, 1 )--Finds all volumes and adds to labels
  
  local  cents = {}
  local nBubbles = labels.Data:max().value
  local volTemp = Double:new()
  if nBubbles>2000 then--Too many bubbles takes too long
    error("Too many bubbles! " .. nBubbles)
  end
  print('# of bubbles: ' .. nBubbles)

  
  for i = 1,nBubbles do
    volTemp = labels:volume(i,i).value
    if volTemp < maxVol and volTemp>0.001 then--Only include bubble if it's volume is between 0.001 and maxVol
      AVS:FIELD_THRESHOLD( labels.Data, temp, i, i, 255, 0)--Put just current bubble into temp
      local clip = chDist(temp,clipDist)--Check if bubbles is too close to convex hull boundary
      if clip then--onyl include if it isn't too close
        local tempCent = field:new()
        AVS:FIELD_CENTER_DOT(temp, tempCent)
        table.insert( cents, {id=i, cent=tempCent})--Put cnetre of bubble into table
      end
    end
    if i == (nBubbles-nBubbles%2)/2 then--For debugging
      print("Halfway through bubble checking")
    end
  end
  
  wm.Scan[7].Adjust = mAdjust--Set up scan to show unflooded bubbles
  wm.Scan[7].Transform = mTrans
  
  for i = 1, #cents do
    AVS:FIELD_THRESHOLD(labels.Data, temp, cents[i].id, cents[i].id, 255, 0)
    wm.Scan[7].Data:add(temp)--Put unflooded bubbles into scan 7
  end
  
  return cents
end

function hasBled(inScan)--Checks to see if the flood fill has bled out by checking total flooded volume
  local inThr = Scan:new()
  local labels = Scan:new()
  inThr.Adjust = mAdjust
  inThr.Transform = mTrans
  labels.Adjust = mAdjust
  labels.Transform = mTrans
  local dummy = field:new()
  local volTemp,volTot = Double:new()
  volTot = 0
  
  AVS:FIELD_THRESHOLD(inScan.Data, inThr.Data, 5000, 5000, 255, 0)--Thresholf off only flooded values
  
  AVS:FIELD_TO_INT(inThr.Data, labels.Data)
  AVS:FIELD_LABEL( labels.Data, labels.Data, dummy, AVS.FIELD_LABEL_3D, 1 )--Find volumes of flooded values
  
  local nBubbles = labels.Data:max().value
  
  for i=1,nBubbles do--add up all the volumes
    volTemp = labels:volume(i,i).value
    volTot = volTot+volTemp
  end

  
  local bleed = false
  if volTot>=bleedVol then--if the volume is greater than bleedVol then return true
    bleed = true
  elseif volTot>maxVolTot then
    maxVolTot=volTot--if it hasn't bled, store the volume (for debugging)
  end
  
  return bleed
end


testpatientpack = getRandScan()--get a random scan
print("Going to use " .. testpatientpack)
loadpack(basefolder .. [[Pack\]] .. testpatientpack)--Load the pack
wm.Scan[1].Description = testpatientpack

--set master adjust and transform.
mAdjust = wm.Scan[1].Adjust
mTrans = wm.Scan[1].Transform

--Make convex hull
local lungCH
lungCH = lungConvexHull()

--threshold into scan 4
wm.Scan[4].Adjust = mAdjust
wm.Scan[4].Transform = mTrans
wm.Scan[5].Adjust = mAdjust
wm.Scan[5].Transform = mTrans
AVS:FIELD_THRESHOLD( wm.Scan[1].Data, wm.Scan[4].Data, 1270, 3000 )

--Smooth and put in scan 5
AVS:FIELD_OPS( wm.Scan[4].Data, wm.Scan[5].Data, 2, AVS.FIELD_OPS_Smooth )--Smoothing filter
local tempScan = Field:new()
tempScan = wm.Scan[5].Data
AVS:FIELD_THRESHOLD( tempScan, wm.Scan[5].Data, 1, 255)--Flatten
AVS:FIELD_OPS(wm.Scan[5].Data, wm.Scan[5].Data, 5, AVS.FIELD_OPS_Closing)--Closing filter to minimise # of volumes

spine = Scan:new()
local delin = wm.Delineation.SC
spine = wm.Scan[1]:burn(delin, 255, true)
AVS:FIELD_OPS(spine.Data, spine.Data, AVS.FIELD_OPS_SignedDist)
local spineBig = Scan:new()
AVS:FIELD_THRESHOLD(spine.Data, spineBig.Data,1)
wm.Scan[5].Data:add(spineBig.Data)

--Mask with CH and put in scan 6. Shrink by 1cm
mask(5,3,6,1)


local cents = bubbleCent(6)--find bubbles and centres
print("Found centres, starting flooding")
local halfway = (#cents-#cents%2)/2
wm.Scan[8] = wm.Scan[1]:copy()
--mask(8,3,8)
for i=1,#cents do--floodfill from centre of each bubble
  local safe = wm.Scan[8]:copy()
  AVS:FLOODFILL(cents[i].cent, wm.Scan[8].Data, 1200, 3, 5000)
  local bleed = hasBled(wm.Scan[8])--check for bleed
  collectgarbage()
  if bleed then--if it has bled, don't include
    wm.Scan[8] = safe:copy()
  end
  if i == halfway then
    print("Halway through flooding")
  end
end
print("Finished flooding")

wm.Scan[1].Description = [[Original CT Scan]]
wm.Scan[2].Description = [[Dose]]
wm.Scan[3].Description = [[Mask]]
wm.Scan[4].Description = [[Thresholded]]
wm.Scan[5].Description = [[Smoothed]]
wm.Scan[6].Description = [[Masked]]
wm.Scan[7].Description = [[Calcifications]]
wm.Scan[8].Description = [[Original with calcifications flooded to 5000]]

print("Finished")