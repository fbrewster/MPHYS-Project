--Check Heart Calcifications
--[[
Frank Brewster
Finds calcifications and then checks what percentage of these are within the heart delineation
]]--


basefolder = [[D:\MPHYS\Data\]]
clipDist = 0.3--acceptabel distance from the centre of a bubble to the CH boundary
maxVol = 5--largest volumes considered a calcification
bleedVol = 10--volume above which a fill is considred to have bled
maxVolTot = 0--to keep track of how big the bleed volumes are
floodVol = 0
lastFloodVol = 0

local log = io.open(basefolder .. [[log.txt]])
local fileList = io.open(basefolder .. [[list.txt]])
local fNames = {}
local index  = 1

if fileList then
  for line in fileList:lines() do
    fNames[index] = unpack(line:split(" "))
    index = index +1
  end
else
  error("No file name list could be found")
end
io.close(fileList)


function TScan:setup()
  self.adjust = mAdjust
  self.Transform = mTrans
end


function lungConvexHull()--Makes a convex hull from the delineation of both lungs
  
  local lungs = wm.Delineation.Both_Lungs or wm.Delineation.Lungs or wm.Delineation.Both_Lung --get lung delineation
  --local lungs = wm.Delineation.body or wm.Delineation.Body or wm.Delineation.Both_Lungs or wm.Delineation.Lungs or wm.Delineation.Both_Lung
   
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



function chDist(bubble, cDist, chIn)--Checks if a bubble is too close to the convex hull in scan 3
  local chTemp = chIn:copy()
  
  AVS:FIELD_OPS(chTemp.Data, chTemp.Data, AVS.FIELD_OPS_SignedDist)--turn ch into a distance measure
  local centre = Field:new()
  AVS:FIELD_CENTER_DOT(bubble, centre)-- put the centre of bubble into centre
  local distF = Field:new()
  AVS:FIELD_SAMPLE(centre, chTemp.Data, distF)-- sample the field ch at all the points in centre and output to dist
  local dist = distF:getvalue(0)--extract distance value
  local distCent = (dist.value - 127)/100--convert to cm from CH boundary
  local distAbs = math.abs(distCent)--take absolute value to include inside and out
  
  local clip = true
  
  if distAbs<cDist then-- if the distance from the centre of the bubble to the ch is less than clipDist the return false
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
    log:write("Too many bubbles! " .. nBubbles .. " Pack: " .. currentpatientpack .. "\n")
    AVS:FIELD_OPS( inSc.Data, inSc.Data, 1, AVS.FIELD_OPS_Smooth )
    AVS:FIELD_TO_INT( inSc.Data, labels.Data)
    AVS:FIELD_LABEL( labels.Data, labels.Data, dummy, AVS.FIELD_LABEL_3D, 1 )
    nBubbles = labels.Data:max().value
    if nBubbles>2000 then
      wm.Scan[4] = inSc:copy()
      error("Tried smoothing again but still too many bubbles. " .. nBubbles)
    else
      log:write("Smoothing again worked" .. "\n")
    end
  end
  log:write('# of bubbles: ' .. nBubbles .. "\n")

  
  for i = 1,nBubbles do
    volTemp = labels:volume(i,i).value
    if volTemp < maxVol and volTemp>0.001 then
      AVS:FIELD_THRESHOLD( labels.Data, temp, i, i, 255, 0)
      local clip = chDist(temp,clipDist,cHull)
      if true then
        local tempCent = field:new()
        AVS:FIELD_CENTER_DOT(temp, tempCent)
        table.insert( cents, {id=i, cent=tempCent})
      end
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

function hasBled2(inScan)
  local vol = inScan:volume(5000,5000).value
  local bled = false
  local newVol = vol-lastFloodVol
  if newVol>maxVol then
    bled = true
  else
    lastFloodVol = vol
  end
  return bled
end


for fileNum=1,index-1 do
  local currentpatientpack = fNames[fileNum]
  log:write(currentpatientpack .. "\n")
  loadpack(basefolder ..[[PackWithHearts\]] .. currentpatientpack)
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
  masked = mask(orign,cHull,0.7)
  local hist = masked:histogram(masked,650,3000,3000,256)
  hist.cumulative = true
  local lowThresh = hist:percentile(99)
  local highThresh = hist:max()

  --threshold into scan 4
  local threshed = Scan:new()
  threshed:setup()
  AVS:FIELD_THRESHOLD( orign.Data, threshed.Data, lowThresh.value, highThresh.value )
  scans[2] = threshed

  --Smooth and put in scan 5
  local smoothed = Scan:new()
  smoothed:setup()
  AVS:FIELD_OPS( threshed.Data, smoothed.Data, 2, AVS.FIELD_OPS_Smooth )--Smoothing filter
  local tempScan = smoothed.Data:copy()
  AVS:FIELD_THRESHOLD( tempScan, smoothed.Data, 1, 255)--Flatten
  AVS:FIELD_OPS(smoothed.Data, smoothed.Data, 5, AVS.FIELD_OPS_Closing)--Closing filter to minimise # of volumes
  print("smoothed")

  --spine burn
  spine = Scan:new()
  spine:setup()
  local delin = wm.Delineation.SC or wm.Delineation.cord or wm.Delineation.SCanal or wm.Delineation.SCANAL or wm.Delineation.SCord or wm.Delineation.sc or wm.Delineation.[S Canal] or wm.Delineation.[s canal]
  spine = orign:burn(delin, 255, true)
  local kern = Field:new("field 3D float", 3,3,3)
  spineWarp = Scan:new()
  spineWarp:setup()

  AVS:MAKE_KERNEL(kern, kern, AVS.MAKE_KERNEL_Gradient, AVS.MAKE_KERNEL_Rectangle, AVS.MAKE_KERNEL_None, 1, 0,0,1)
  local spineBig = Scan:new()
  spineBig:setup()
  AVS:FIELD_OPS(spine.Data, spine.Data, 15, AVS.FIELD_OPS_Smooth3)
  AVS:FIELD_THRESHOLD(spine.Data, spineBig.Data,1)
  smoothed.Data:add(spineBig.Data)
  scans[3] = smoothed

  --Mask with CH and put in scan 6. Shrink by 1cm
  local maskSmooth = Scan:new()
  maskSmooth:setup()
  maskSmooth = mask(smoothed,cHull,0.7)
  scans[4] = maskSmooth


  local cents = field:new()
  local bubTot = Scan:new()
  bubTot:setup()
  bubs = {}
  local shcHull = cHull:copy()
  AVS:FIELD_OPS(shcHull.Data, shcHull.Data, AVS.FIELD_OPS_SignedDist)--make singed distance image
  tempScan = shcHull.Data:copy()
  AVS:FIELD_THRESHOLD(tempScan, shcHull.Data, 197)--delete everything shrink cm in
  scans[1] = shcHull
  cents,bubTot,bubs = bubbleCent(maskSmooth,shcHull:copy())--find bubbles and centres
  scans[5] = bubTot
  local flood = masked:copy()
  local allMask = Scan:new()
  allMask:setup()

  for i=1,#cents do--floodfill from centre of each bubble
    local currentBub = bubs[i]:copy()
    local bubMask = Scan:new()
    bubMask:setup()
    AVS:FIELD_OPS(currentBub.Data, currentBub.Data, 4, AVS.FIELD_OPS_Smooth3)
    AVS:FIELD_MASK(orign.Data, currentBub.Data, bubMask.Data)
    local bubHist = bubMask:histogram(bubMask,650,3000,3000,256)
    bubHist.cumulative = true
    local floodThresh = bubHist:percentile(80).value
    if floodThresh<1100 then 
      floodThresh = 1135 
    end
    allMask:add(bubMask)
    
    local safe = flood:copy()
    AVS:FLOODFILL(cents[i].cent, flood.Data, floodThresh, 3, 5000)
    
    local bleed = hasBled2(flood)--check for bleed
    collectgarbage()
    if bleed then--if it has bled, don't include
      flood = safe:copy()
      floodThresh = bubHist:percentile(83).value
      AVS:FLOODFILL(cents[i].cent, flood.Data, floodThresh, 3, 5000)
      local bleed2 = hasBled2(flood)
      if bleed2 then
        flood = safe:copy()
      end
    end
  end

  local null = field:new()
  local bubFlood = Scan:new()
  bubFlood:setup()
  local totBubFlood = Scan:new()
  totBubFlood:setup()
  AVS:FIELD_THRESHOLD( flood.Data, totBubFlood.Data, 5000,5000)
  AVS:FIELD_OPS(totBubFlood.Data, totBubFlood.Data, 1, AVS.FIELD_OPS_Smooth3)
  tempScan = totBubFlood.Data:copy()
  AVS:FIELD_THRESHOLD(tempScan, totBubFlood.Data, 1, 255)
  totBubFlood:add(spineBig)
  AVS:FIELD_TO_INT( totBubFlood.Data, bubFlood.Data )
  AVS:FIELD_LABEL( bubFlood.Data, bubFlood.Data, null, AVS.FIELD_LABEL_3D, 1)
  local doseScan = wm.Scan[2]:copy()
  local thisBub = Scan:new()
  thisBub:setup()
  local finalBubs = Scan:new()
  finalBubs:setup()

  for i=1,#cents do
    AVS:FIELD_THRESHOLD(bubFlood.Data, thisBub.Data, i, i)
    if thisBub.Data:max().value>0 and thisBub:volume().value<10 then
      local bubHist = thisBub:histogram(doseScan, 255, 255, 100, 100)
      bubHist.cumulative = true
      local median = bubHist:percentile(50).value
      if median>0 then
        table.insert(medianDose, median)
      end
      finalBubs:add(thisBub)
    end
    collectgarbage()
    thisBub.Data:clear()
  end

  collectgarbage()
end



io.close(log)