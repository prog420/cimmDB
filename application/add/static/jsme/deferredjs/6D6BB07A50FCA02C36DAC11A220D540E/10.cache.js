$wnd.jsme.runAsyncCallback10('var e8="data-selenium-id";function f8(){this.lb=ox("file");this.lb[gi]="gwt-FileUpload";this.a=new g8;this.a.c=this;if(-1==this.hb){var a=this.lb,b=4096|(this.lb.__eventBits||0);vD();eE(a,b)}else this.hb|=4096}x(399,380,Gm,f8);_.Ae=function(a){var b;a:{b=this.a;switch(tD(a.type)){case 1024:if(!b.a){b.b=!0;b=!1;break a}break;case 4096:if(b.b){b.a=!0;var c=b.c.lb,d=rx(di,!0);c.dispatchEvent(d);b.a=!1;b.b=!1}}b=!0}b&&xE(this,a)};_.a=null;x(400,1,{});function g8(){}x(401,400,{},g8);_.a=!1;_.b=!1;\n_.c=null;function h8(a){var b=$doc.createElement(yi);or(il,b.tagName);this.lb=b;this.b=new K0(this.lb);this.lb[gi]="gwt-HTML";J0(this.b,a,!0);S0(this)}x(405,406,Gm,h8);function i8(a,b){var c,d;c=$doc.createElement(Kl);d=$doc.createElement(wl);d[Fh]=a.a.a;d.style[Sl]=a.b.a;var e=(oD(),pD(d));c.appendChild(e);nD(a.d,c);JE(a,b,d)}function j8(){HF.call(this);this.a=(KF(),RF);this.b=(SF(),VF);this.e[bi]=gd;this.e[ai]=gd}x(454,396,Bm,j8);\n_.Ve=function(a){var b;b=qx(a.lb);(a=NE(this,a))&&this.d.removeChild(qx(b));return a};function k8(a){try{a.s=!1;var b,c,d;d=a.db;c=a.Y;d||(a.lb.style[Tl]=fj,a.Y=!1,a.gf());b=a.lb;b.style[qj]=0+(Ky(),Fk);b.style[Fl]=md;N2(a,sZ(yx($doc)+(xx()-lx(a.lb,mk)>>1),0),sZ(zx($doc)+(wx()-lx(a.lb,lk)>>1),0));d||((a.Y=c)?(a.lb.style[ji]=Mk,a.lb.style[Tl]=Ul,at(a.cb,200)):a.lb.style[Tl]=Ul)}finally{a.s=!0}}function l8(a,b){var c;c=(new D1(a)).rd.hg();c.lb.setAttribute(e8,"jsa_clipboard/button/"+b);return c}\nfunction m8(a){var b;b=l8("Close (ESC)","close");tE(b,new n8(a),(Pz(),Pz(),Qz));return b}\nfunction o8(){A2();var a,b,c,d,e;Z2.call(this,(r3(),s3),null,!0);this.kj();this._=!0;this.lb.setAttribute(e8,"jsa_clipboard/window");this.R=!0;a=new h8(this.e);this.d=new cH;this.d.lb.setAttribute(e8,"jsa_clipboard/text_area");mE(this.d,od);iE(this.d,od);r2(this,"400px");e=new j8;e.lb.style[dj]=od;e.e[bi]=10;c=(KF(),LF);e.a=c;i8(e,a);i8(e,this.d);this.c=new ZF;this.c.e[bi]=20;for(b=this.ij(),c=0,d=b.length;c<d;++c)a=b[c],WF(this.c,a);i8(e,this.c);F2(this,e);P2(this,!1);tE(this.d,new p8(this),(DA(),\nDA(),EA));this.jj()}x(803,804,AZ,o8);_.ij=function(){return z(oH,o,51,[m8(this)])};_.jj=function(){var a=this.d;a.lb.readOnly=!0;var b=nE(a.lb)+"-readonly";hE(a.Ie(),b,!0)};_.kj=function(){q3(this.E.b,"Copy")};_.gf=function(){Y2(this);this.lb.style[Yl]=ud};_.c=null;_.d=null;_.e="Press Ctrl-C (Command-C on Mac) or right click (Option-click on Mac) on the selected text to copy it, then paste into another program.";function p8(a){this.a=a}x(806,1,{},p8);\n_.pe=function(a){27==(a.a.keyCode||0)&&H2(this.a,!1)};_.a=null;function n8(a){this.a=a}x(807,1,{},n8);_.ge=function(){H2(this.a,!1)};_.a=null;function q8(a){this.a=a}x(808,1,{},q8);_.Od=function(){oE(this.a.d.lb,!0);CF(this.a.d,!0);var a=this.a.d,b;b=mx(a.lb,Rl).length;if(0<b&&a.gb){if(0>b)throw new dT("Length must be a positive integer. Length: "+b);if(b>mx(a.lb,Rl).length)throw new dT("From Index: 0  To Index: "+b+"  Text Length: "+mx(a.lb,Rl).length);try{a.lb.setSelectionRange(0,0+b)}catch(c){}}};\n_.a=null;function r8(a){var b;b=l8(a.a,"accept");tE(b,new s8(a),(Pz(),Pz(),Qz));return b}function t8(a){a.e="Paste the text to import into the text area below.";a.a="Accept";q3(a.E.b,"Paste")}function u8(a){A2();o8.call(this);this.b=a}x(810,803,AZ,u8);_.ij=function(){return z(oH,o,51,[r8(this),m8(this)])};_.jj=function(){iE(this.d,"150px")};_.kj=function(){t8(this)};_.gf=function(){Y2(this);this.lb.style[Yl]=ud;ax((Yw(),Zw),new v8(this))};_.a=null;_.b=null;function w8(a){A2();u8.call(this,a)}\nx(809,810,AZ,w8);_.ij=function(){var a;return z(oH,o,51,[r8(this),(a=new f8,a.lb.setAttribute(e8,"jsa_clipboard/button/browse_upload"),tE(a,new x8(this),(F_(),F_(),G_)),a),m8(this)])};_.jj=function(){iE(this.d,"150px");AL(new y8(this),this.d)};_.kj=function(){t8(this);this.e+=" Or drag and drop a file on it."};function x8(a){this.a=a}x(811,1,{},x8);_.fe=function(a){var b,c;b=new FileReader;a=(c=a.a.target,c.files[0]);z8(b,new A8(this));b.readAsText(a)};_.a=null;function A8(a){this.a=a}\nx(812,1,{},A8);_.vg=function(a){$G(this.a.a.d,a)};_.a=null;function y8(a){this.a=a;this.b=new B8(this);this.c=this.d=1}x(813,566,{},y8);_.a=null;function B8(a){this.a=a}x(814,1,{},B8);_.vg=function(a){this.a.a.d.lb[Rl]=null!=a?a:m};_.a=null;function s8(a){this.a=a}x(818,1,{},s8);_.ge=function(){if(this.a.b){var a=this.a.b,b;b=new rK(a.a,0,mx(this.a.d.lb,Rl));FQ(a.a.a,b.a)}H2(this.a,!1)};_.a=null;function v8(a){this.a=a}x(819,1,{},v8);_.Od=function(){oE(this.a.d.lb,!0);CF(this.a.d,!0)};_.a=null;\nx(820,1,jn);_.Zd=function(){var a,b;a=new C8(this.a);void 0!=$wnd.FileReader?b=new w8(a):b=new u8(a);t2(b);k8(b)};function C8(a){this.a=a}x(821,1,{},C8);_.a=null;x(822,1,jn);_.Zd=function(){var a;a=new o8;var b=this.a,c,d;$G(a.d,b);c=(d=BT(b,"\\r\\n|\\r|\\n|\\n\\r"),d.length);1>=c&&(c=~~(b.length/16));iE(a.d,20*(10>c+1?c+1:10)+Fk);ax((Yw(),Zw),new q8(a));t2(a);k8(a)};function z8(a,b){a.onload=function(a){b.vg(a.target.result)}}Z(803);Z(810);Z(809);Z(821);Z(806);Z(807);Z(808);Z(818);Z(819);Z(811);Z(812);\nZ(813);Z(814);Z(405);Z(454);Z(399);Z(400);Z(401);R(uZ)(10);\n//@ sourceURL=10.js\n')
