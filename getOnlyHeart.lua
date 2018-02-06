--Extract only calcifications in the heart.

basefolder = [[D:\MPHYS\Data\]]
dataFileName = [[PackWithHearts\]]
outFolder = [[Calcifications\]]
xdrFolder = outFolder .. [[heartMasks\]]
origxdrFolder = outFolder .. [[masks\]]

local function scandir(directory)
  local i, t, popen = 0, {}, io.popen;
  for filename in popen('dir "'..directory..'" /o:n /b /a-d'):lines() do
    i = i + 1;
    t[i] = filename;
  end
  return t;
end

function TScan:setup()--sets adjust and transform
  self.adjust = mAdjust;
  self.Transform = mTrans;
end

local function contains (tab, val)
    for index, value in ipairs(tab) do
        if value == val then
            return true;
        end
    end

    return false;
end

local function median( t )
  local temp={}

  -- deep copy table so that when we sort it, the original is unchanged
  -- also weed out any non numbers
  for k,v in pairs(t) do
    if type(v) == 'number' then
      table.insert( temp, v )
    end
  end

  table.sort( temp )

  -- If we have an even number of table elements or odd.
  if math.fmod(#temp,2) == 0 then
    -- return mean value of middle two elements
    return ( temp[#temp/2] + temp[(#temp/2)+1] ) / 2
  else
    -- return middle element
    return temp[math.ceil(#temp/2)]
  end
end

--Get list of scans with heart delin
--log = io.open(basefolder .. outFolder .. [[log.txt]], 'a');--log file used for error handling and debugging
volFile = io.open(basefolder .. outFolder .. [[heartVols.csv]], 'w');
volFile:write( "ID,mVol \n" );

local fxdrNames = scandir(basefolder .. origxdrFolder);
local fHeartNames = scandir(basefolder .. dataFileName);
local fNames = {};

for i=1,#fHeartNames do
  fHeartNames[i]=fHeartNames[i]:gsub(".pack", ".xdr");
  if contains(fxdrNames,fHeartNames[i]) then
    table.insert( fNames, { id = i, pack = fHeartNames[i]:gsub(".xdr",".pack") } );
  end
end
collectgarbage()

for i=1,#fNames do
  currentpatientpack = fNames[i].pack;
--load xdr and heart delin
  loadpack( basefolder .. dataFileName .. currentpatientpack );
  mAdjust = wm.Scan[1].Adjust;
  mTrans = wm.Scan[1].Transform;
  local origBubs = field:new();
  AVS:READ_XDR( origBubs, basefolder .. origxdrFolder .. currentpatientpack:gsub( ".pack", ".xdr" ) );

--mask off heart
  local heartDelin = wm.Delineation.Heart or wm.Delineation.heart;
  local heartBurn = Scan:new();
  local masked = Scan:new();
  heartBurn:setup();
  masked:setup();
  heartBurn = wm.Scan[1]:burn( heartDelin, 255, 0 );
  AVS:FIELD_MASK( origBubs, heartBurn.data, masked.Data );
  

--discard small vols
  local nOfBubs = masked.Data:max().value;
  local thisBub = Scan:new();
  local finalBubs = Scan:new();
  thisBub:setup();
  finalBubs:setup();
  local val = 1;
  local vols = {};
  
  for i=1,nOfBubs do
    AVS:FIELD_THRESHOLD( masked.Data, thisBub.Data, i, i, val )--take one volume and give it unique pixel value
    local thisVol = thisBub:volume().value;
    if thisBub.Data:max().value>0 and thisVol>0.01 then--if the bubble exists and isn't too small
      finalBubs:add( thisBub );--add it to the final bubbles
      table.insert( vols , thisVol );
      val = val + 1;--move to next pixel value id
    end
    thisBub.Data:clear();
    collectgarbage();
  end

--output mVols
  local mVol = 0;
  if #vols>1 then
    mVol = median( vols );
  elseif #vols==1 then
    mVol = vols[1];
  end
  volFile:write( currentpatientpack:gsub( ".pack", "" ) .. "," .. mVol .. " \n" );
  volFile:flush();
  
  collectgarbage();
end

io.close(volFile);