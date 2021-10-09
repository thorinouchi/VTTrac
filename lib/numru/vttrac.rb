# VTTrac, Velocimetry by Template Tracking, Ruby interface
#
# author: Takeshi Horinouchi
# date:  2021-05-26
#
# Copyright 2021 (C) Takeshi Horinouchi. All rights reserved.
# 
# LICENSE: see ../../LICENSE.txt
# 
# This is a C wrapper for the C-library VTTrac by the author.
# A typical usage is like this:
#
#    # initialize a VTTrac object:
#    vtt = VTTrac.new(z, t: t, tunit: "s", zmiss: -999.0)
#
#    # set up tracking parameters:
#    vtt.setup( nsx: 20, nsy: 20, ... )  # ... : further keyword-option setting
#
#    # conduct tracking:
#    count, x, y, vx, vy, score, zss, score_ary =
#                vtt.trac(tid, x, y, out_subimage: true, out_score_ary: true)


require "narray"
require "narray_miss"
require "numru/vttrac_core"

module NumRu
  class VTTrac < VTTracCore

    require "numru/vttrac_version"  # required here after inheriting VTTracCore

    SCORE_XCOR = 1
    SCORE_NCOV = 2

    attr_reader :nt, :setup, :vxhw, :vyhw
    
    # Object initialization (VTTrac)
    # 
    # You need to call #setup subsequently before calling #trac
    # * z [Array or 2D NArrays or 3D NArray]: the "images" used for tracking
    #   (when 3D, last dimension should be time)
    # * t [1D NArray of nil]: time for each image (nil->0,1,2,..)
    # * tunit [String or nil]: like "second". Can be nil (non-dimensional)
    # * zmiss [Float] : missing value in z (nil if there is no missing)
    def initialize(z, t: nil, tunit: nil, zmiss: nil)
      case z
      when Array
        @nt = z.length
        nxy = z[0].length
        (1...@nt).each{|it|
          ze = z[it]
          if ze.length!=nxy
            raise(ArgumentError, "Shapes of the arrays are not uniform")
          end
          if !ze.is_a?(NArray)
            raise(ArgumentError,"z element must be a NArray (now #{ze.class})")
          end
        }
      when NArray
        @nt = z.shape[-1]
      else
        raise(ArgumentError,"Unsupported class for z [#{z.class}]")
      end

      case t
      when nil
        t = NArray.double(@nt).indgen!
      when NArray
        raise(ArgumentError,"invalid time axis length") if t.length != @nt
      else
        raise(ArgumentError,"Unsupported class for time [#{z.class}]")
      end

      case tunit
      when nil
        tunit = ""
      end

      if zmiss
        chk_zmiss = true
      else
        chk_zmiss = false
        zmiss = 0.0
      end
      
      super(z, chk_zmiss, zmiss, t, tunit)
      @setup = false
    end

    # Setup for tracking
    # 
    # * nsx, nsy (MANDADORY, Integer): submimage x & y sizes (x:1st, y:2nd dim)
    # * vxhw, vyhw (Float): (either v[xy]hw or i[xy]hw are MANDATORY)
    #   search velocity range half sizes to set i[xy]hw.
    #   Seach at least to cover +-v?hw around the first guess or previous step.
    #   (the result can be outside the range)
    # * ixhw, iyhw (Integer): (either v[xy]hw or i[xy]hw are MANDATORY)
    #   max displacement fro template match (can be set indirecly thru v[xy]hw)
    # * subgrid [true/false, default:true]: whether to conduct subgrid tracking
    # * subgrid_grid [true/false, default:true]: Whether subgrid peak finding is by gaussian
    # * itstep [Integer>0, default: 1]: step of t's used (skip if >1)
    # * ntrack [Integer>0, default: 2]: max tracking times from initial loc
    # * score_method {String]: "xcor" for cross-correlation (so far,
    #   this is the only method)
    # * score_th0 [Float]: the minimum score required for the 1st tracking
    # * score_th1 [Float]: the minimum score required for subsequent tracking
    # * vxch [nil or Float]: if non-nil, the max tolerant vx change between
    #   two consecutive traking.
    # * vych [nil or Float]: if non-nil, the max tolerant vy change between
    #   two consecutive traking.
    def setup(nsx:, nsy:, # mandatory
              vxhw: nil, vyhw: nil, ixhw: nil, iyhw: nil, # semi-mandatory
              subgrid: true, subgrid_gaus: false, itstep: 1, ntrac: 2,
              score_method: :xcor,  score_th0: 0.8,  score_th1: 0.7,
              vxch: nil, vych: nil, peak_inside: true, peak_inside_th: 0.03,
              min_contrast: nil, use_init_temp: false)
      if vxhw
        raise("v[xy]hw and i[xy]hw must not be set simultaneously") if ixhw
        raise("vxhw and vyhw must be set simultaneously") if !vyhw
        set_ixyhw_from_v(vxhw, vyhw)
      elsif ixhw
        raise("ixhw and iyhw must be set simultaneously") if !iyhw
        set_ixyhw_directly(ixhw, iyhw)
      else
        raise("either i[xy]hw or v[xy]hw must be specified")
      end
      @vxhw = vxhw
      @vyhw = vyhw
      @ixhw = ixhw
      @iyhw = iyhw
      vxch = vxch || -999.0    # <=0 for nil (not to set)
      vych = vych || -999.0    # <=0 for nil (not to set)

      case score_method
      when :xcor
        score_method = SCORE_XCOR
      when :ncov
        score_method = SCORE_NCOV
      else
        raise("unsupported method #{score_method}")
      end

      if !peak_inside
        peak_inside_th = -1.0  # negative, meaning unused in the C library
      end
      if !min_contrast
        min_contrast = -1.0   # negative, meaning unused in the C library
      end

      set_basic(nsx, nsy,
                itstep, ntrac)
      set_optional(subgrid, subgrid_gaus, score_method,
               score_th0, score_th1,
               peak_inside_th, min_contrast,
               vxch, vych, use_init_temp)

      @setup = true
      nil
    end

    # Conduct tracking
    #
    # * tid [Integer or integer NArray]: tracking initial t indices.
    #   If Numeric, treated as a constant common for all of x,y,etc.
    # * x [float/integer NArray]: Tracking initial template-center x location (index-based; non-integer for subgrid)
    # * y [float/integer NArray]: Tracking initial template-center y location (index-based; non-integer for subgrid)
    # * vxg [optional float NArray]: vx initial guess (if nil, set to 0)
    # * vyg [optional float NArray]: vy initial guess (if nil, set to 0)
    # * out_subimage: [false: default, or true] if true, returens the
    #   subimages along the track zss (see below)
    # 
    # REMARK
    # The shapes of tid, x, y, vxg, vyg must be the same (rank can be > 1).
    # This shape shall be expressed as *sh in what follows.
    # e.g. if if the shape of the input arrays are [k,l], the shapes of
    # the outputs shall be like, count: [k,l], tid: [ntrac+1,k,l]
    # 
    # RETURN VALUE
    # * An Array containing the following
    #   * count (int NArray[*sh])  The number of successful tracking for each initial template.
    #   * tid (int NArray[ntrac+1,*sh])  time indices of the trajectories (starting from the initial ones)
    #   * x (float NArray[ntrac+1,*sh])  x locations of the trajectories (starting from the initial ones)
    #   * y (float NArray[ntrac+1,*sh])  y locations of the trajectories (starting from the initial ones)
    #   * vx (float NArray[ntrac,*sh])  Derived x-velocity
    #   * vy (float NArray[ntrac,*sh])  Derived y-velocity
    #   * score (float NArray[ntrac,*sh])  Scores along the trajectory (max values, possibly at subgrid)
    #   * zss  (nil if !out_subimage, or 4D NArray[nsx,nsy, ntrac+1,*sh])
    #     the entire subimages along the tracks  (including the initial one)
    #   * score_ary (nil if !out_subimage, or 4D NArray[(x-sliding size),(y-sliding size), ntrac+1,*sh])
    #     the entire scores
    # 
    def trac(tid, x, y, vxg: nil, vyg: nil,
             out_subimage: false, out_score_ary: false)
      raise(ArgumentError,"Need to call #setup in advance") unless @setup
      sh = x.shape
      if tid.is_a?(Integer)
        tid = NArray.int(*sh).fill!(tid)
      else
        raise(ArgumentError,"Shape miss-match (x)") if tid.shape != sh
      end
      raise(ArgumentError,"Shape miss-match (y)") if y.shape != sh 
      if vxg.nil?
        vxg = NArray.float(*sh)
      else
        raise(ArgumentError,"Shape miss-match (vxg)") if vxg.shape != sh
      end
      if vyg.nil?
        vyg = NArray.float(*sh)
      else
        raise(ArgumentError,"Shape miss-match (vyg)") if vyg.shape != sh
      end
      if !tid.is_a?(NArray) || tid.typecode!=NArray::INT
        raise("tid != int NArray")
      end

      result = super(tid, x, y, vxg, vyg, out_subimage, out_score_ary)
      result[2..6] = result[2..6].map{|na| NArrayMiss.to_nam(na, na.ne(rmiss))}
      if sh.length >= 2
        # reshape outputs based on the shape of inputs
        result[0] = result[0].reshape(*sh)
        result[1..8] = result[1..8].map{|na|
          if na
            newshape = na.shape[0..-2] + sh
            na.reshape(*newshape)
          else
            na
          end
        }
      end
      result
    end
    
  end
end

#################################################
#### test part
#################################################

if $0 == __FILE__
  include NumRu
  include NMath
  iws = ( ARGV[0] || 1 ).to_i
  ifl = 1
  nt = 30
  nx = 100
  ny = 100
  t = NArray.float(nt).indgen!
  tunit = ""
  zmiss = -999.0
  xg = NArray.sfloat(nx).indgen!.newdim(-1,-1)
  yg = NArray.sfloat(ny).indgen!.newdim(0,-1)
  tg = t.newdim(0,0)
  k = 2*PI / 10.0
  #k = 2*PI / 9.3
  c = 1.2
  z = sin( k*(xg-c*tg) )*cos( k*(yg-c*tg) )
  vtt = VTTrac.new(z, t: t, tunit: tunit, zmiss: zmiss)
  ntrac = nt-1
  nsx = 5
  nsy = 5
  #nsx = 8
  #nsy = 8
  vtt.setup(nsx: nsx, nsy: nsy, vxhw: 1.8, vyhw: 1.8, ntrac: ntrac,
            score_method: :xcor, peak_inside: false)
  n = 6
  tid0 = NArray.int(n)
  x0 = NArray.float(n).indgen!*2.5 + 7.5
  y0 = NArray.float(n).indgen!*1 + 10.5
  count, tid, x, y, vx, vy, score, zss, score_ary =
               vtt.trac(tid0, x0, y0, out_subimage: true, out_score_ary: true)
#exit

  p "** count", count
  p "** x, y", x, y
  p "** vx, vy", vx, vy
  p "** vx, vy  mean", vx.mean(0), vy.mean(0)
  p "** score", score
  #p "** zss", zss[2,2,false]
  #p "** score_ary max", score_ary.max(0,1)
  p "** score_ary", score_ary
  p "///** score_ary", score_ary[true,true,8,0]
#exit
  require "numru/dcl"
  DCL.swpset("iwidth",800)
  DCL.swpset("iheight",800)
  DCL.swpset("lalt",true)
  DCL.swpset("lwait",false)
  DCL.swpset("ifl",ifl)
  DCL.gropn(iws)
  DCL.sgscmn(5)
  DCL.uepset("ltone",true)
  DCL.uepset("icolor1",40)

  (0..ntrac).each do |it|
    DCL.grfrm
    DCL.grsvpt(0.15, 0.85, 0.15, 0.85)
    DCL.grswnd(0, nx-1, 0, ny-1)
    #DCL.grswnd(5, 21, 5, 21)
    DCL.grstrf
    DCL.usdaxs
    DCL.uegtla(-1,1,0.1)
    DCL.uetonc(z[true,true,it])
    #DCL.uetonc(z[5..21,5..21,it])
    (0...x.shape[-1]).each{|k|
      xc = x[it,k]
      yc = y[it,k]
      x1 = xc - nsx/2 - 0.5
      x2 = x1 + nsx
      y1 = yc - nsy/2 - 0.5
      y2 = y1 + nsy
      bx = [x1,x1,x2,x2,x1]
      by = [y1,y2,y2,y1,y1]
      DCL.sgplzu(bx, by, 1, 1 + 20*(k%2) )
    }
    sleep 1.0 if it==1
    sleep 0.3 if it>1
  end

  vxg = x0.dup.fill(5.0)  # half wavelength right
  vyg = x0.dup.fill(5.0)  # half wavelength up     --> there are same features
  count, tid, x, y, vx, vy, score, zss, score_ary =
                vtt.trac(tid0, x0, y0, out_subimage: true, out_score_ary: true,
                vxg: vxg, vyg: vyg  )
  (0..10).each do |it|
    DCL.grfrm
    DCL.grsvpt(0.15, 0.85, 0.15, 0.85)
    DCL.grswnd(0, nx-1, 0, ny-1)
    #DCL.grswnd(5, 21, 5, 21)
    DCL.grstrf
    DCL.usdaxs
    DCL.uxsttl("t","if vxg&vyg are set to too far",0.0)
    DCL.uegtla(-1,1,0.1)
    DCL.uetonc(z[true,true,it])
    #DCL.uetonc(z[5..21,5..21,it])
    (0...x.shape[-1]).each{|k|
      xc = x[it,k]
      yc = y[it,k]
      x1 = xc - nsx/2 - 0.5
      x2 = x1 + nsx
      y1 = yc - nsy/2 - 0.5
      y2 = y1 + nsy
      bx = [x1,x1,x2,x2,x1]
      by = [y1,y2,y2,y1,y1]
      DCL.sgplzu(bx, by, 1, 1 + 20*(k%2) )
    }
    sleep 0.3
  end

  vtt.use_init_temp = true
  count, tid, x, y, vx, vy, =
               vtt.trac(tid0, x0, y0)

  (0..ntrac).each do |it|
    DCL.grfrm
    DCL.grsvpt(0.15, 0.85, 0.15, 0.85)
    DCL.grswnd(0, nx-1, 0, ny-1)
    #DCL.grswnd(5, 21, 5, 21)
    DCL.grstrf
    DCL.usdaxs
    DCL.uxsttl("t","when use_init_temp is true",0.0)
    DCL.uegtla(-1,1,0.1)
    DCL.uetonc(z[true,true,it])
    #DCL.uetonc(z[5..21,5..21,it])
    (0...x.shape[-1]).each{|k|
      xc = x[it,k]
      yc = y[it,k]
      x1 = xc - nsx/2 - 0.5
      x2 = x1 + nsx
      y1 = yc - nsy/2 - 0.5
      y2 = y1 + nsy
      bx = [x1,x1,x2,x2,x1]
      by = [y1,y2,y2,y1,y1]
      DCL.sgplzu(bx, by, 1, 1 + 20*(k%2) )
    }
    sleep 1.0 if it<=1
    sleep 0.3 if it>1
    #gets
  end
  p "** vx, vy  mean", vx.mean(0), vy.mean(0)

  DCL.grcls
end
