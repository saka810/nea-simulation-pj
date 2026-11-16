// NEA GeoSim アイコン変換・比較一覧生成。作業者: Codex
// 実行: node export_icons.cjs [sharpパッケージの絶対パス]
const fs = require('node:fs/promises');
const path = require('node:path');
const crypto = require('node:crypto');
const SIZES = [16,20,24,32,40,48,64,128,256];
const escape = s => s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
async function main(sharp) {
 const root = path.resolve(__dirname,'..');
 const designs = JSON.parse(await fs.readFile(path.join(__dirname,'デザイン.json'),'utf8'));
 for(const d of ['PNG','ICO']) await fs.mkdir(path.join(root,d),{recursive:true});
 // 元絵を切り抜かず、縦横比を保って正方形へ収める。白背景と余白を保持。
 const resize = (source,n) => sharp(source).resize(n,n,{fit:'contain',background:'#FFFCF8',kernel:'lanczos3'}).png().toBuffer();
 const fontfile = 'C:/Windows/Fonts/YuGothM.ttc';
 const overlays = [];
 async function label(text,x,y,size=20,color='#102C54') {
  const b = await sharp({text:{text:'<span foreground="'+color+'">'+escape(text)+'</span>',font:'Yu Gothic '+size,fontfile,dpi:72,rgba:true}}).png().toBuffer();
  overlays.push({input:b,left:x,top:y});
 }
 const manifest=[];
 await label('NEA GeoSim（ネアシム）｜線画アイコン 8案',40,26,32);
 await label('各案：大きな見本 ＋ 16 / 32 / 48px の実寸比較（明・暗背景）',40,78,19);
 await label('白い背景を保持。暗色欄も同じPNGを配置。画面の100%表示で縮小欄を確認。',40,108,16);
 for(let i=0;i<designs.length;i++){
  const d=designs[i], file=path.join(root,'原寸PNG',d.id+'.png');
  const original=await fs.readFile(file), meta=await sharp(original).metadata();
  await fs.writeFile(path.join(root,'PNG',d.id+'.png'),await resize(original,512));
  const images=await Promise.all(SIZES.map(n=>resize(original,n)));
  // ICOディレクトリと各PNGフレームを明示的に組み立てる（9サイズを確実に格納）。
  const head=Buffer.alloc(6+16*SIZES.length); head.writeUInt16LE(1,2);head.writeUInt16LE(SIZES.length,4);
  let offset=head.length;
  for(let j=0;j<SIZES.length;j++){
   const n=SIZES[j],at=6+j*16;
   head[at]=n===256?0:n;head[at+1]=n===256?0:n;
   head.writeUInt16LE(1,at+4);head.writeUInt16LE(32,at+6);
   head.writeUInt32LE(images[j].length,at+8);head.writeUInt32LE(offset,at+12);
   offset+=images[j].length;
  }
  const ico=Buffer.concat([head,...images]);await fs.writeFile(path.join(root,'ICO',d.id+'.ico'),ico);
  // 書いたICOの全フレームを読み返し、サイズとPNGの復号を検査。
  if(ico.readUInt16LE(4)!==9)throw Error('ICOのサイズ数が不正');
  for(let j=0;j<9;j++){
   const at=6+j*16,sz=ico.readUInt32LE(at+8),start=ico.readUInt32LE(at+12);
   const m=await sharp(ico.subarray(start,start+sz)).metadata();
   if(m.width!==SIZES[j]||m.height!==SIZES[j])throw Error('ICOフレームの寸法が不正');
  }
  const x=30+(i%4)*390,y=160+Math.floor(i/4)*590;
  overlays.push({input:Buffer.from('<svg width="374" height="566"><rect x="0" y="0" width="374" height="566" rx="12" fill="white" stroke="#DCE3E8"/></svg>'),left:x,top:y});
  await label(d.id+'  '+d.name,x+16,y+18,21);
  overlays.push({input:await resize(original,276),left:x+49,top:y+62});
  await label('明背景',x+16,y+357,15); for(const [j,n] of [16,32,48].entries()) await label(n+'px',x+[95,173,255][j],y+357,15);
  const positions=[x+113,x+191,x+273];
  for(let j=0;j<3;j++){const n=[16,32,48][j];overlays.push({input:await resize(original,n),left:positions[j]-Math.floor(n/2),top:y+398-Math.floor(n/2)});}
  overlays.push({input:Buffer.from('<svg width="342" height="96"><rect width="342" height="96" rx="6" fill="#202A38"/></svg>'),left:x+16,top:y+437});
  await label('暗背景',x+28,y+470,14,'#FFFFFF');
  for(let j=0;j<3;j++){const n=[16,32,48][j];overlays.push({input:await resize(original,n),left:positions[j]-Math.floor(n/2),top:y+485-Math.floor(n/2)});}
  manifest.push({id:d.id,name:d.name,original:{width:meta.width,height:meta.height,sha256:crypto.createHash('sha256').update(original).digest('hex')},png:{width:512,height:512},ico_sizes:SIZES,verified:true});
 }
 await sharp({create:{width:1620,height:1370,channels:3,background:'#F1F4F5'}}).composite(overlays).png().toFile(path.join(root,'一覧.png'));
 await fs.writeFile(path.join(__dirname,'検証結果.json'),JSON.stringify(manifest,null,2),'utf8');
 return manifest;
}
module.exports={main};
if(require.main===module){main(require(process.argv[2]||'sharp')).then(x=>console.log(x.length+'案を出力・検証しました')).catch(e=>{console.error(e);process.exitCode=1;});}
