// Generated by CoffeeScript 1.4.0
(function() {
  var Animation, BrowserSupport, Gravity, Spring, Spring2, Tween, TweenForce, TweenGravity, TweenLinear, TweenSpring, TweenSpring2,
    __bind = function(fn, me){ return function(){ return fn.apply(me, arguments); }; },
    __hasProp = {}.hasOwnProperty,
    __extends = function(child, parent) { for (var key in parent) { if (__hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; };

  Tween = (function() {

    Tween.properties = {};

    function Tween(options) {
      this.options = options != null ? options : {};
      this.next = __bind(this.next, this);

      this.init = __bind(this.init, this);

    }

    Tween.prototype.init = function() {
      return this.t = 0;
    };

    Tween.prototype.next = function(step) {
      if (this.t > 1) {
        this.t = 1;
      }
      this.currentT = this.t;
      return this.t += step;
    };

    return Tween;

  })();

  TweenLinear = (function(_super) {

    __extends(TweenLinear, _super);

    function TweenLinear() {
      this.next = __bind(this.next, this);

      this.init = __bind(this.init, this);
      return TweenLinear.__super__.constructor.apply(this, arguments);
    }

    TweenLinear.prototype.init = function() {
      return TweenLinear.__super__.init.apply(this, arguments);
    };

    TweenLinear.prototype.next = function(step) {
      var t;
      TweenLinear.__super__.next.call(this, step);
      t = this.currentT;
      return [t, t];
    };

    return TweenLinear;

  })(Tween);

  TweenForce = (function(_super) {

    __extends(TweenForce, _super);

    function TweenForce() {
      this.next = __bind(this.next, this);

      this.init = __bind(this.init, this);
      return TweenForce.__super__.constructor.apply(this, arguments);
    }

    TweenForce.prototype.init = function() {
      return TweenForce.__super__.init.apply(this, arguments);
    };

    TweenForce.prototype.next = function(step) {
      var gravity, t, v;
      TweenForce.__super__.next.call(this, step);
      t = this.currentT;
      gravity = 1;
      v = 2 * t - gravity * t * t;
      return [t, v];
    };

    return TweenForce;

  })(Tween);

  TweenGravity = (function(_super) {

    __extends(TweenGravity, _super);

    function TweenGravity() {
      this.next = __bind(this.next, this);

      this.curve = __bind(this.curve, this);

      this.init = __bind(this.init, this);

      this.length = __bind(this.length, this);

      this.gravityValue = __bind(this.gravityValue, this);

      this.bounceValue = __bind(this.bounceValue, this);

      this.duration = __bind(this.duration, this);
      return TweenGravity.__super__.constructor.apply(this, arguments);
    }

    TweenGravity.properties = {
      bounce: {
        min: 0,
        max: 99,
        "default": 40
      },
      gravity: {
        min: 1,
        max: 4000,
        "default": 1000
      },
      duration: {
        editable: false
      }
    };

    TweenGravity.prototype.duration = function() {
      return Math.round(1000 * 1000 / this.options.gravity * this.length());
    };

    TweenGravity.prototype.bounceValue = function() {
      var bounce;
      bounce = this.options.bounce / 100;
      bounce = Math.min(bounce, 99);
      return bounce;
    };

    TweenGravity.prototype.gravityValue = function() {
      return this.options.gravity / 100;
    };

    TweenGravity.prototype.length = function() {
      var L, b, bounce, curve, gravity;
      bounce = this.bounceValue();
      gravity = this.gravityValue();
      b = Math.sqrt(2 / gravity);
      curve = {
        a: -b,
        b: b,
        H: 1
      };
      while (curve.H > 0.001) {
        L = curve.b - curve.a;
        curve = {
          a: curve.b,
          b: curve.b + L * bounce,
          H: curve.H * bounce * bounce
        };
      }
      return curve.b;
    };

    TweenGravity.prototype.init = function() {
      var L, b, bounce, curve, gravity, _results;
      TweenGravity.__super__.init.apply(this, arguments);
      L = this.length();
      gravity = this.gravityValue();
      gravity = gravity * L * L;
      bounce = this.bounceValue();
      b = Math.sqrt(2 / gravity);
      this.curves = [];
      curve = {
        a: -b,
        b: b,
        H: 1
      };
      this.curves.push(curve);
      _results = [];
      while (curve.b < 1 && curve.H > 0.001) {
        L = curve.b - curve.a;
        curve = {
          a: curve.b,
          b: curve.b + L * bounce,
          H: curve.H * bounce * bounce
        };
        _results.push(this.curves.push(curve));
      }
      return _results;
    };

    TweenGravity.prototype.curve = function(a, b, H, t) {
      var L, t2;
      L = b - a;
      t2 = (2 / L) * t - 1 - (a * 2 / L);
      return t2 * t2 * H - H + 1;
    };

    TweenGravity.prototype.next = function(step) {
      var bounce, curve, gravity, i, t, v;
      TweenGravity.__super__.next.call(this, step);
      t = this.currentT;
      bounce = this.options.bounce / 100;
      gravity = this.options.gravity;
      i = 0;
      curve = this.curves[i];
      while (!(t >= curve.a && t <= curve.b)) {
        i += 1;
        curve = this.curves[i];
        if (!curve) {
          break;
        }
      }
      if (!curve) {
        v = 1;
      } else {
        v = this.curve(curve.a, curve.b, curve.H, t);
      }
      return [t, v];
    };

    return TweenGravity;

  })(Tween);

  TweenSpring = (function(_super) {

    __extends(TweenSpring, _super);

    function TweenSpring() {
      this.next = __bind(this.next, this);
      return TweenSpring.__super__.constructor.apply(this, arguments);
    }

    TweenSpring.properties = {
      frequency: {
        min: 0,
        max: 100,
        "default": 15
      },
      friction: {
        min: 1,
        max: 1000,
        "default": 100
      },
      anticipationStrength: {
        min: 0,
        max: 1000,
        "default": 115
      },
      anticipationSize: {
        min: 0,
        max: 99,
        "default": 10
      },
      duration: {
        min: 100,
        max: 4000,
        "default": 1000
      }
    };

    TweenSpring.prototype.next = function(step) {
      var A, At, a, angle, b, decal, frequency, friction, frictionT, s, t, v, y0, yS,
        _this = this;
      TweenSpring.__super__.next.call(this, step);
      t = this.currentT;
      frequency = Math.max(1, this.options.frequency);
      friction = Math.pow(20, this.options.friction / 100);
      s = this.options.anticipationSize / 100;
      decal = Math.max(0, s);
      frictionT = (t / (1 - s)) - (s / (1 - s));
      if (t < s) {
        A = function(t) {
          var M, a, b, x0, x1;
          M = 0.8;
          x0 = s / (1 - s);
          x1 = 0;
          b = (x0 - (M * x1)) / (x0 - x1);
          a = (M - b) / x0;
          return (a * t * _this.options.anticipationStrength / 100) + b;
        };
        yS = (s / (1 - s)) - (s / (1 - s));
        y0 = (0 / (1 - s)) - (s / (1 - s));
        b = Math.acos(1 / A(yS));
        a = (Math.acos(1 / A(y0)) - b) / (frequency * (-s));
      } else {
        A = function(t) {
          return Math.pow(friction / 10, -t) * (1 - t);
        };
        b = 0;
        a = 1;
      }
      At = A(frictionT);
      angle = frequency * (t - s) * a + b;
      v = 1 - (At * Math.cos(angle));
      return [t, v, At, frictionT, angle];
    };

    return TweenSpring;

  })(Tween);

  TweenSpring2 = (function(_super) {

    __extends(TweenSpring2, _super);

    function TweenSpring2() {
      this.next = __bind(this.next, this);
      return TweenSpring2.__super__.constructor.apply(this, arguments);
    }

    TweenSpring2.properties = {
      frequency: {
        min: 0,
        max: 100,
        "default": 15
      },
      friction: {
        min: 1,
        max: 1000,
        "default": 100
      },
      duration: {
        min: 100,
        max: 4000,
        "default": 1000
      }
    };

    TweenSpring2.prototype.next = function(step) {
      var A, At, At2, Ax, angle, frequency, friction, t, v,
        _this = this;
      TweenSpring2.__super__.next.call(this, step);
      t = this.currentT;
      frequency = Math.max(1, this.options.frequency);
      friction = Math.pow(20, this.options.friction / 100);
      A = function(t) {
        return 1 - Math.pow(friction / 10, -t) * (1 - t);
      };
      At = A(t);
      At2 = A(1 - t);
      Ax = (Math.cos(t * 2 * 3.14 - 3.14) / 2) + 0.5;
      Ax = Math.pow(Ax, this.options.friction / 100);
      angle = frequency * t;
      v = Math.cos(angle) * Ax;
      return [t, v, Ax, -Ax];
    };

    return TweenSpring2;

  })(Tween);

  BrowserSupport = (function() {

    function BrowserSupport() {}

    BrowserSupport.transform = function() {
      return this.withPrefix("transform");
    };

    BrowserSupport.keyframes = function() {
      if (document.body.style.webkitAnimation !== void 0) {
        return "-webkit-keyframes";
      }
      if (document.body.style.mozAnimation !== void 0) {
        return "-moz-keyframes";
      }
      return "keyframes";
    };

    BrowserSupport.withPrefix = function(property) {
      var prefix;
      prefix = this.prefixFor(property);
      if (prefix !== '') {
        return "-" + (prefix.toLowerCase()) + "-" + property;
      }
      return property;
    };

    BrowserSupport.prefixFor = function(property) {
      var k, prefix, prop, propArray, propertyName, _i, _j, _len, _len1, _ref;
      propArray = property.split('-');
      propertyName = "";
      for (_i = 0, _len = propArray.length; _i < _len; _i++) {
        prop = propArray[_i];
        propertyName += prop.substring(0, 1).toUpperCase() + prop.substring(1);
      }
      _ref = ["Webkit", "Moz"];
      for (_j = 0, _len1 = _ref.length; _j < _len1; _j++) {
        prefix = _ref[_j];
        k = prefix + propertyName;
        if (document.body.style[k] !== void 0) {
          return prefix;
        }
      }
      return '';
    };

    return BrowserSupport;

  })();

  Animation = (function() {

    Animation.index = 0;

    Animation.prototype.tweenClass = "TweenLinear";

    function Animation(el, frames, options) {
      var _base;
      this.el = el;
      this.frames = frames != null ? frames : {};
      this.options = options != null ? options : {};
      this._keyframes = __bind(this._keyframes, this);

      this.start = __bind(this.start, this);

      this.tween = __bind(this.tween, this);

      (_base = this.options).duration || (_base.duration = 1000);
    }

    Animation.prototype.tween = function() {
      this._tween || (this._tween = eval("new " + this.tweenClass + "(this.options)"));
      return this._tween;
    };

    Animation.prototype.start = function() {
      var animation, k, keyframes, name, prefix, property, propertyName, style, v, _results;
      name = "anim_" + Animation.index;
      Animation.index += 1;
      keyframes = this._keyframes(name);
      style = document.createElement('style');
      style.innerHTML = keyframes;
      document.head.appendChild(style);
      animation = {
        name: name,
        duration: this.options.duration + 'ms',
        timingFunction: 'linear',
        fillMode: 'forwards'
      };
      _results = [];
      for (k in animation) {
        v = animation[k];
        property = "animation-" + k;
        prefix = BrowserSupport.prefixFor(property);
        propertyName = prefix + "Animation" + k.substring(0, 1).toUpperCase() + k.substring(1);
        _results.push(this.el.style[propertyName] = v);
      }
      return _results;
    };

    Animation.prototype._keyframes = function(name) {
      var args, css, dValue, frame0, frame1, isTransform, k, newValue, oldValue, properties, step, t, transform, unit, v, value;
      this.tween().init();
      step = 0.01;
      frame0 = this.frames[0];
      frame1 = this.frames[100];
      css = "@" + (BrowserSupport.keyframes()) + " " + name + " {\n";
      while (args = this.tween().next(step)) {
        t = args[0], v = args[1];
        transform = '';
        properties = {};
        for (k in frame1) {
          value = frame1[k];
          value = parseFloat(value);
          oldValue = frame0[k] || 0;
          dValue = value - oldValue;
          newValue = oldValue + (dValue * v);
          unit = '';
          isTransform = false;
          if (k === 'translateX' || k === 'translateY' || k === 'translateZ') {
            unit = 'px';
            isTransform = true;
          } else if (k === 'rotateX' || k === 'rotateY' || k === 'rotateZ') {
            unit = 'deg';
            isTransform = true;
          } else if (k === 'scaleX' || k === 'scaleY' || k === 'scale') {
            isTransform = true;
          }
          if (isTransform) {
            transform += "" + k + "(" + newValue + unit + ") ";
          } else {
            properties[k] = newValue;
          }
        }
        css += "" + (t * 100) + "% {\n";
        if (transform) {
          css += "" + (BrowserSupport.transform()) + ": " + transform + ";\n";
        }
        for (k in properties) {
          v = properties[k];
          css += "" + k + ": " + v + ";\n";
        }
        css += " }\n";
        if (t >= 1) {
          break;
        }
      }
      css += "}\n";
      return css;
    };

    return Animation;

  })();

  Spring = (function(_super) {

    __extends(Spring, _super);

    Spring.prototype.tweenClass = "TweenSpring";

    Spring.properties = TweenSpring.properties;

    function Spring(el, from, to, options) {
      this.el = el;
      this.from = from;
      this.to = to;
      this.options = options != null ? options : {};
      this.frames = {
        0: this.from,
        100: this.to
      };
      Spring.__super__.constructor.call(this, this.el, this.frames, this.options);
    }

    return Spring;

  })(Animation);

  Spring2 = (function(_super) {

    __extends(Spring2, _super);

    Spring2.prototype.tweenClass = "TweenSpring2";

    Spring2.properties = TweenSpring2.properties;

    function Spring2(el, from, to, options) {
      this.el = el;
      this.from = from;
      this.to = to;
      this.options = options != null ? options : {};
      this.frames = {
        0: this.from,
        100: this.to
      };
      Spring2.__super__.constructor.call(this, this.el, this.frames, this.options);
    }

    return Spring2;

  })(Animation);

  Gravity = (function(_super) {

    __extends(Gravity, _super);

    Gravity.prototype.tweenClass = "TweenGravity";

    Gravity.properties = TweenGravity.properties;

    function Gravity(el, from, to, options) {
      this.el = el;
      this.from = from;
      this.to = to;
      this.options = options != null ? options : {};
      this.frames = {
        0: this.from,
        100: this.to
      };
      this.options.duration = this.tween().duration();
      Gravity.__super__.constructor.call(this, this.el, this.frames, this.options);
    }

    return Gravity;

  })(Animation);

  this.Dynamics = {
    Spring: Spring,
    Spring2: Spring2,
    Gravity: Gravity,
    BrowserSupport: BrowserSupport
  };

}).call(this);
