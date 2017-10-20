--Quantify Calcifications
--[[
Frank Brewster
Creates convex hull around both lungs delineation and then thresholds to find calcifications
]]--

local testpatientpack = [[201103328.pack]]
basefolder = [[C:\Users\Frank\MPHYS\Data\]]
clipDist = 0.1

function getRandScan()--Gets a random scan name from the list in ..\Data
  local file = io.open(basefolder..[[list.txt]])
  local fNames = {}
  local i = 1
  local fNamesSplit = {}

  if file then
    for line in file:lines() do
        fNames[i] = unpack(line:split(" "))
        i=i+1
    end
  end
  
  local nOfFiles = i-1  
  local index = math.random(nOfFiles)  
  local randScan = fNames[index]
  
  return randScan
end

function lungConvexHull()--Makes a convex hull from the delineation of both lungs
  
  local lungs = wm.Delineation.Both_Lungs or wm.Delineation.Lungs or wm.Delineation.Both_Lung --get lung delineation
  --[[if (lungs==nil) then --find lung delineation if it has a weird name
    local answer
    io.write("Error: No lung delineation found. Is it named something else (y or n)? ")
    --io.flush
    answer = io.read(1)
    if (answer=='y') then
      local lungName
      io.write("Please enter it: ")
      --io.flush
      lungName=io.read(15)
      lungs=wm.Delineation(lungName)
    else
      error("Cannot proceed without a lung delineation")
    end
  end]]--
    
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
  
  wm.Scan[3] = wm.Scan[1]:burn(newDel, 255, false)
  wm.Scan[3].Description = "Lung Convex Hull Delineation"
  
  return newDel
end


function mask(input, cHull, output, shrink)--masks scan[input] with scan[cHull], shrunk by "shrink" and outputs to scan[output]
  shrink = -1*shrink
  --Exclude data outside of CH
  wm.Scan[output].Adjust = wm.Scan[input].Adjust
  wm.Scan[output].Transform = wm.Scan[input].Transform
  --need to shrink scan
  AVS:FIELD_OPS(wm.Scan[cHull].Data, wm.Scan[cHull].Data, AVS.FIELD_OPS_SignedDist)--make singed distance image
  AVS:FIELD_THRESHOLD(wm.Scan[cHull].Data, wm.Scan[cHull].Data, shrink)--delete everything 3cm in
  AVS:FIELD_MASK(wm.Scan[input].Data, wm.Scan[cHull].Data, wm.Scan[output].Data)--delete outside of new shrunk mask
  wm.Scan[cHull].Description = "Convex hull shrunk by 3cm"
  wm.Scan[output].Description = "Masked and thresholded data"
end


function excludeVols(inSc, outSc, maxVol)
  wm.Scan[outSc].Adjust = wm.Scan[inSc].Adjust
  wm.Scan[outSc].Transform = wm.Scan[inSc].Transform
  local dummy = field:new()
  local labels = Scan:new()
  labels.Adjust = wm.Scan[inSc].Adjust
  labels.Transform = wm.Scan[inSc].Transform
  local temp = Field:new()
  local clipDist = 0.1
  
  AVS:FIELD_TO_INT( wm.Scan[inSc].Data, labels.Data)
  AVS:FIELD_LABEL( labels.Data, labels.Data, dummy, AVS.FIELD_LABEL_3D, 1 )
  
  local vols = {}
  local nBubbles = labels.Data:max().value
  if nBubbles>2000 then
    error("Too many bubbles! " .. nBubbles)
  end
  print('# of bubbles: ' .. nBubbles)
  local volTemp
  
  for i = 1,nBubbles do
    volTemp = labels:volume(i,i).value
    if volTemp < maxVol and volTemp>0.01 then
      --AVS:FIELD_THRESHOLD( labels.Data, temp, i, i, 255, 0)
      --local clip = chDist(temp,clipDist)
      --print(clip)
      --if clip then
        table.insert( vols, { id=i, vol=volTemp } )
      --end
      --end
    end
  end
  
  
  
  local function compareVols(a,b)
    return a.vol > b.vol
  end
  table.sort(vols, compareVols)
  
  

  for j=1, #vols do
    AVS:FIELD_THRESHOLD( labels.Data, temp, vols[j].id, vols[j].id, 255, 0)
    if j==1 then
      --local centre = temp:center()
      --print("Largest bubble located at " .. centre)
      print(temp:center())
    end
    wm.Scan[outSc].Data:add(temp)
  end
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
  chOut.Data:distxfm()
  ch.Data:clear()
  ch.Data = chIn.Data
  chOut.Data:tobyte()
  ch.Data:add(chOut.Data)
  local centre = Field:new()
  AVS:FIELD_CENTER_DOT(bubble, centre)-- put the centre of bubble into centre
  local dist = Field:new()
  AVS:FIELD_SAMPLE(centre, ch.Data, dist)-- sample the field ch at all the points in centre and output to dist
    print(dist:getvalue(1,1,1))
  
  if distAbs<=clipDist then-- if the distance from the centre of the bubble to the ch is less than clipDist the return false
    local clip = false
  else
    local clip = true
  end
  
  return clip
end
testpatientpack = getRandScan()
print("Going to use " .. testpatientpack)
loadpack(basefolder .. [[Pack\]] .. testpatientpack)
wm.Scan[1].Description = testpatientpack

--Make convex hull
local lungCH
lungCH = lungConvexHull()

--threshold into scan 4
--local toThresh = Scan:new()
--toThresh.Data = wm.Scan[1].Data
wm.Scan[4].Adjust = wm.Scan[1].Adjust
wm.Scan[4].Transform = wm.Scan[1].Transform
wm.Scan[5].Adjust = wm.Scan[1].Adjust
wm.Scan[5].Transform = wm.Scan[1].Transform
--AVS:FIELD_OPS( toThresh.Data, toThresh.Data, 5, AVS.FIELD_OPS_Smooth)
AVS:FIELD_THRESHOLD( wm.Scan[1].Data, wm.Scan[4].Data, 1135, 1650 )
AVS:FIELD_OPS( wm.Scan[4].Data, wm.Scan[5].Data, 1, AVS.FIELD_OPS_Smooth )
local tempScan = Field:new()
--tempScan.Adjust = wm.Scan[4].Adjust
--tempScan.Transform = wm.Scan[4].Transform
tempScan = wm.Scan[5].Data
AVS:FIELD_THRESHOLD( tempScan, wm.Scan[5].Data, 1, 255)
wm.Scan[4].Description = "Thresholded Scan"

--[[--Exclude data outside of CH
wm.Scan[5] = wm.Scan[1]
--need to shrink scan
AVS:FIELD_OPS(wm.Scan[3].Data, wm.Scan[3].Data, AVS.FIELD_OPS_SignedDist)--make singed distance image
AVS:FIELD_THRESHOLD(wm.Scan[3].Data, wm.Scan[3].Data, -3)--delete everything 3cm in
AVS:FIELD_MASK(wm.Scan[4].Data, wm.Scan[3].Data, wm.Scan[5].Data)--delete outside of new shrunk mask
wm.Scan[3].Description = "Convex hull shrunk by 3cm"
wm.Scan[5].Description = "Masked and thresholded data"]]--

AVS:FIELD_OPS(wm.Scan[5].Data, wm.Scan[5].Data, 5, AVS.FIELD_OPS_Closing)
mask(5,3,6,2.8)--mask scan 4 

excludeVols(6,7,3)

--print("Pack: " .. testpatientpack)
print("Finished")