package org.reactome.immport.ws.test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.httpclient.HttpClient;
import org.apache.commons.httpclient.HttpMethod;
import org.apache.commons.httpclient.HttpStatus;
import org.apache.commons.httpclient.methods.GetMethod;
import org.apache.commons.httpclient.methods.PostMethod;
import org.apache.commons.httpclient.methods.RequestEntity;
import org.apache.commons.httpclient.methods.StringRequestEntity;
import org.junit.Test;
import org.reactome.immport.ws.model.requests.GSMForVOs;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;

public class WSTests {
    protected final String HOST_URL = "http://localhost:8076/immportws";
    protected final String HTTP_POST = "Post";
    protected final String HTTP_GET = "Get";

    
    @Test
    public void testGetStudy() throws Exception {
        String url = HOST_URL + "/study/SDY212";
        testGetMethod(url);
    }
    
    private void testGetMethod(String url) throws Exception {
        System.out.println(url);
        String rtn = callHttp(url, HTTP_GET, null);
        outputJSON(rtn);
    }
    
    @Test
    public void testRepository() throws Exception {
        String url = HOST_URL + "/repository/GEO";
        testGetMethod(url);
    }
    
    @Test
    public void testQueryStudiesForVO() throws Exception {
    	String url = HOST_URL + "/study/vaccine/VO_0000642";
    	String rtn = callHttp(url, HTTP_GET, null);
    	outputJSON(rtn);
    }
    
    @Test
    public void testQueryStudiesForVOs() throws Exception {
    	String url = HOST_URL + "/study/vaccine";
    	String query = "VO_0000044,VO_0000642,VO_0004809,VO_0000045,VO_0000046,VO_0000047";
    	String rtn = callHttp(url, HTTP_POST, query);
    	outputJSON(rtn);
    }
    
    @Test
    public void testQueryTimesForVO() throws Exception {
    	String url  = HOST_URL + "/collectionTimes/vaccine/VO_0000642";
    	String rtn = callHttp(url,HTTP_GET, null);
    	outputJSON(rtn);
    }
    
    @Test
    public void testQueryTimesForVOs() throws Exception {
    	String url  = HOST_URL + "/collectionTimes/vaccine";
    	String query = "VO_0004809,VO_0000047,VO_0000642,VO_0000044,VO_0000045"+"\n"+"male,female";
    	String rtn = callHttp(url, HTTP_POST, query);
    	outputJSON(rtn);
    }
    
    @Test
    public void testAnalyzeBiosamples() throws Exception {
    	String url = HOST_URL + "/analysis/geneExpression";
    	String rtn = callHttp(url,HTTP_POST, "{\n"
    			+ "  \"GSMids\": [\"GSM2712782\",\"GSM2712783\",\"GSM2712784\",\"GSM2712785\",\"GSM2712786\",\"GSM2712787\",\"GSM2712788\",\"GSM2712789\",\"GSM2712790\",\"GSM2712791\",\"GSM2712792\",\"GSM2712793\",\"GSM2712794\",\"GSM2712795\",\"GSM2712796\",\"GSM2712797\",\"GSM2712798\",\"GSM2712799\",\"GSM2712800\",\"GSM2712801\",\"GSM2712802\",\"GSM2712803\",\"GSM2712804\",\"GSM2712805\",\"GSM2712806\",\"GSM2712807\",\"GSM2712808\",\"GSM2712809\",\"GSM2712810\",\"GSM2712811\",\"GSM2712812\",\"GSM2712813\",\"GSM2712814\",\"GSM2712815\",\"GSM2712816\",\"GSM2712817\",\"GSM2712818\",\"GSM2712819\",\"GSM2712820\",\"GSM2712821\",\"GSM2712822\",\"GSM2712823\",\"GSM2712824\",\"GSM2712825\",\"GSM2712826\",\"GSM2712827\",\"GSM2712828\",\"GSM2712829\",\"GSM2712830\",\"GSM2712831\",\"GSM2712832\",\"GSM2712833\",\"GSM2712834\",\"GSM2712835\",\"GSM2712836\",\"GSM2712837\",\"GSM2712838\",\"GSM2712839\",\"GSM2712840\",\"GSM2712841\",\"GSM2712842\",\"GSM2712843\",\"GSM2712844\",\"GSM2712845\",\"GSM2712846\",\"GSM2712847\",\"GSM2712848\",\"GSM2712849\",\"GSM2712850\",\"GSM2712851\",\"GSM2712852\",\"GSM2712853\",\"GSM2712854\",\"GSM2712855\",\"GSM2712856\",\"GSM2712857\",\"GSM2712858\",\"GSM2712859\",\"GSM2712860\",\"GSM2712861\",\"GSM2712862\",\"GSM2712863\",\"GSM2712864\",\"GSM2712865\",\"GSM2712866\",\"GSM2712867\",\"GSM2712868\",\"GSM2712869\",\"GSM2712870\",\"GSM2712871\",\"GSM2712872\",\"GSM2712873\",\"GSM2712874\",\"GSM2712875\",\"GSM2712876\",\"GSM2712877\",\"GSM2712878\",\"GSM2712879\",\"GSM1441830\",\"GSM1441831\",\"GSM1441832\",\"GSM1441833\",\"GSM1441834\",\"GSM1441835\",\"GSM1441836\",\"GSM1441837\",\"GSM1441838\",\"GSM1441839\",\"GSM1441840\",\"GSM1441841\",\"GSM1441842\",\"GSM1441843\",\"GSM1441844\",\"GSM1441845\",\"GSM1441846\",\"GSM1441847\",\"GSM1441848\",\"GSM1441849\",\"GSM1441850\",\"GSM1441851\",\"GSM1441852\",\"GSM1441853\",\"GSM1441854\",\"GSM1441855\",\"GSM1441856\",\"GSM1441857\",\"GSM1441858\",\"GSM1441859\",\"GSM1441860\",\"GSM1441861\",\"GSM1441862\",\"GSM1441863\",\"GSM1441864\",\"GSM1441865\",\"GSM1441866\",\"GSM1441867\",\"GSM1441868\",\"GSM1441869\",\"GSM1441870\",\"GSM1441871\",\"GSM1441872\",\"GSM1441873\",\"GSM1441874\",\"GSM1441875\",\"GSM1441876\",\"GSM1441877\",\"GSM1441878\",\"GSM1441879\",\"GSM1441880\",\"GSM1441881\",\"GSM1441882\",\"GSM1441883\",\"GSM1441884\",\"GSM1441885\",\"GSM1441886\",\"GSM1441887\",\"GSM1441888\",\"GSM1441889\",\"GSM1441890\",\"GSM1441891\",\"GSM1441892\",\"GSM1441893\",\"GSM1441894\",\"GSM1441895\",\"GSM1441896\",\"GSM1441897\",\"GSM1441898\",\"GSM1441899\",\"GSM1441900\",\"GSM1441901\",\"GSM1441902\",\"GSM1441903\",\"GSM1441904\",\"GSM1441905\",\"GSM1441906\",\"GSM1441907\",\"GSM1441908\",\"GSM1441909\",\"GSM1441910\",\"GSM1441911\",\"GSM1441912\",\"GSM1441913\",\"GSM1441914\",\"GSM1441915\",\"GSM1441916\",\"GSM1441917\",\"GSM1441918\",\"GSM1441919\",\"GSM1441920\",\"GSM1441921\",\"GSM1441922\",\"GSM1441923\",\"GSM1441924\",\"GSM1441925\",\"GSM1441926\",\"GSM1441927\",\"GSM1441928\",\"GSM1441929\",\"GSM1441930\",\"GSM1441931\",\"GSM1441932\",\"GSM1441933\",\"GSM1441934\",\"GSM1441935\",\"GSM1441936\",\"GSM1441937\",\"GSM1441938\",\"GSM1441939\",\"GSM1441940\",\"GSM1441941\",\"GSM1441942\",\"GSM1441943\",\"GSM1441944\",\"GSM1441945\",\"GSM1441946\",\"GSM1441947\",\"GSM1441948\",\"GSM1441949\",\"GSM1441950\",\"GSM1441951\",\"GSM1441952\",\"GSM1441953\",\"GSM1441954\",\"GSM1441955\",\"GSM1441956\",\"GSM1441957\",\"GSM1441958\",\"GSM1441959\",\"GSM1441960\",\"GSM1441961\",\"GSM1441962\",\"GSM1441963\",\"GSM1441964\",\"GSM1441965\",\"GSM1441966\",\"GSM1441967\",\"GSM1441968\",\"GSM1441969\",\"GSM1441970\",\"GSM1441971\",\"GSM1441972\",\"GSM1441973\",\"GSM1441974\",\"GSM1441975\",\"GSM1441976\",\"GSM1441977\",\"GSM1441978\",\"GSM1441979\",\"GSM1441980\",\"GSM1441981\",\"GSM1441982\",\"GSM1441983\",\"GSM1441984\",\"GSM1441985\",\"GSM1445822\",\"GSM1445823\",\"GSM1445824\",\"GSM1445825\",\"GSM1445826\",\"GSM1445827\",\"GSM1445828\",\"GSM1445829\",\"GSM1445830\",\"GSM1445831\",\"GSM1445832\",\"GSM1445833\",\"GSM1445834\",\"GSM1445835\",\"GSM1445836\",\"GSM1445837\",\"GSM1445838\",\"GSM1445839\",\"GSM1445840\",\"GSM1445841\",\"GSM1445842\",\"GSM1445843\",\"GSM1445844\",\"GSM1445845\",\"GSM1445846\",\"GSM1445847\",\"GSM1445848\",\"GSM1445849\",\"GSM1445850\",\"GSM1445851\",\"GSM1445852\",\"GSM1445854\",\"GSM1445855\",\"GSM1445857\",\"GSM1445860\",\"GSM1445862\",\"GSM1445864\",\"GSM1445866\",\"GSM1445868\",\"GSM1445869\",\"GSM1445870\",\"GSM1445871\",\"GSM1445872\",\"GSM1445873\",\"GSM1445874\",\"GSM1445875\",\"GSM1445876\",\"GSM1445877\",\"GSM1445878\",\"GSM1445879\",\"GSM1445880\",\"GSM1445881\",\"GSM1445882\",\"GSM1445883\",\"GSM1445884\",\"GSM1445885\",\"GSM1445886\",\"GSM1445887\",\"GSM1445888\",\"GSM1445889\",\"GSM1445890\",\"GSM1445891\",\"GSM1445892\",\"GSM1445893\",\"GSM1445894\",\"GSM1445895\",\"GSM1445896\",\"GSM1445897\",\"GSM1445898\",\"GSM1445899\",\"GSM1445900\",\"GSM1445901\",\"GSM1445902\",\"GSM1445903\",\"GSM1445904\",\"GSM1445905\",\"GSM1445906\",\"GSM1445907\",\"GSM1445908\",\"GSM1445909\",\"GSM1445910\",\"GSM1445911\",\"GSM1445912\",\"GSM1445913\",\"GSM1445914\",\"GSM1445915\",\"GSM1445916\",\"GSM1445917\",\"GSM1445918\",\"GSM1445919\",\"GSM1445920\",\"GSM1445921\",\"GSM1445922\",\"GSM1445923\",\"GSM1445924\",\"GSM1445925\",\"GSM1445926\",\"GSM1445927\",\"GSM1445928\",\"GSM1445929\",\"GSM1445930\",\"GSM1445931\",\"GSM1445932\",\"GSM1445933\",\"GSM1445934\",\"GSM1445935\",\"GSM1445936\",\"GSM1445937\",\"GSM1445938\",\"GSM1445939\",\"GSM1445940\",\"GSM1445941\",\"GSM1445942\",\"GSM1445943\",\"GSM1445944\",\"GSM1445945\",\"GSM1445946\",\"GSM1445947\",\"GSM1445948\",\"GSM1445949\",\"GSM1441474\",\"GSM1441475\",\"GSM1441476\",\"GSM1441477\",\"GSM1441478\",\"GSM1441479\",\"GSM1441480\",\"GSM1441481\",\"GSM1441482\",\"GSM1441483\",\"GSM1441484\",\"GSM1441485\",\"GSM1441486\",\"GSM1441487\",\"GSM1441488\",\"GSM1441489\",\"GSM1441490\",\"GSM1441491\",\"GSM1441492\",\"GSM1441493\",\"GSM1441494\",\"GSM1441495\",\"GSM1441496\",\"GSM1441497\",\"GSM1441498\",\"GSM1441499\",\"GSM1441500\",\"GSM1441501\",\"GSM1441502\",\"GSM1441503\",\"GSM1441504\",\"GSM1441505\",\"GSM1441506\",\"GSM1441507\",\"GSM1441508\",\"GSM1441509\",\"GSM1441510\",\"GSM1441511\",\"GSM1441512\",\"GSM1441513\",\"GSM1441514\",\"GSM1441515\",\"GSM1441516\",\"GSM1441517\",\"GSM1441518\",\"GSM1441519\",\"GSM1441520\",\"GSM1441521\",\"GSM1441522\",\"GSM1441523\",\"GSM1441524\",\"GSM1441525\",\"GSM1441526\",\"GSM1441527\",\"GSM1441528\",\"GSM1441529\",\"GSM1441530\",\"GSM1441531\",\"GSM1441532\",\"GSM1441533\",\"GSM1441534\",\"GSM1441535\",\"GSM1441536\",\"GSM1441537\",\"GSM1441538\",\"GSM1441539\",\"GSM1441540\",\"GSM1441541\",\"GSM1441542\",\"GSM1441543\",\"GSM1441544\",\"GSM1441545\",\"GSM1441546\",\"GSM1441547\",\"GSM1441548\",\"GSM1441549\",\"GSM1441550\",\"GSM1441551\",\"GSM1441552\",\"GSM1441553\",\"GSM1441554\",\"GSM1441555\",\"GSM1441556\",\"GSM1441557\",\"GSM1441558\",\"GSM1441559\",\"GSM1441560\",\"GSM1441561\",\"GSM1441562\",\"GSM1441563\",\"GSM1441564\",\"GSM1441565\",\"GSM1441566\",\"GSM1441567\",\"GSM1441568\",\"GSM1441569\",\"GSM1441570\",\"GSM1441571\",\"GSM1441572\",\"GSM1441573\",\"GSM1441574\",\"GSM1441575\",\"GSM1441576\",\"GSM1441577\",\"GSM1441578\",\"GSM1441579\",\"GSM1441580\",\"GSM1441581\",\"GSM1441582\",\"GSM1441583\",\"GSM1441584\",\"GSM1441585\",\"GSM1441586\",\"GSM1441587\",\"GSM1441588\",\"GSM1441589\",\"GSM1441590\",\"GSM1441591\",\"GSM1441592\",\"GSM1441593\",\"GSM1441594\",\"GSM1441595\",\"GSM1441596\",\"GSM1441597\",\"GSM1441598\",\"GSM1441599\",\"GSM1441600\",\"GSM1441601\",\"GSM1441602\",\"GSM1441603\",\"GSM1441604\",\"GSM1441605\",\"GSM1441606\",\"GSM1441607\",\"GSM1441608\",\"GSM1441609\",\"GSM1441610\",\"GSM1441611\",\"GSM1441612\",\"GSM1441613\",\"GSM1441614\",\"GSM1441615\",\"GSM1441616\",\"GSM1441617\",\"GSM1441618\",\"GSM1441619\",\"GSM1441620\",\"GSM1441621\",\"GSM1441622\",\"GSM1441623\",\"GSM1441624\",\"GSM1441625\",\"GSM1441626\",\"GSM1441627\",\"GSM1441628\",\"GSM1441629\",\"GSM1441630\",\"GSM1441631\",\"GSM1441632\",\"GSM1441633\",\"GSM1441634\",\"GSM1441635\",\"GSM1441636\",\"GSM1441637\",\"GSM1441638\",\"GSM1441639\",\"GSM1441640\",\"GSM1441641\",\"GSM1441642\",\"GSM1441643\",\"GSM1441644\",\"GSM1441645\",\"GSM1441646\",\"GSM1441647\",\"GSM1441648\",\"GSM1441649\",\"GSM1441650\",\"GSM1441651\",\"GSM1441652\",\"GSM1441653\",\"GSM1441654\",\"GSM1441655\",\"GSM1441656\",\"GSM1441657\",\"GSM1441658\",\"GSM1441659\",\"GSM1441660\",\"GSM1441661\",\"GSM1441662\",\"GSM1441663\",\"GSM1441664\",\"GSM1441665\",\"GSM1441666\",\"GSM1441667\",\"GSM1441668\",\"GSM1441669\",\"GSM1441670\",\"GSM1441671\",\"GSM1441672\",\"GSM1441673\",\"GSM1441674\",\"GSM1441675\",\"GSM1441676\",\"GSM1441677\",\"GSM1441678\",\"GSM1441679\",\"GSM1441680\",\"GSM1441681\",\"GSM1441682\",\"GSM1441683\",\"GSM1441684\",\"GSM1441685\",\"GSM1441686\",\"GSM1441687\",\"GSM1441688\",\"GSM1441689\",\"GSM1441690\",\"GSM1441691\",\"GSM1441692\",\"GSM1441693\",\"GSM1441694\",\"GSM1441695\",\"GSM1441696\",\"GSM1441697\",\"GSM1441698\",\"GSM1441699\",\"GSM1441700\",\"GSM1441701\",\"GSM1441702\",\"GSM1441703\",\"GSM1441704\",\"GSM1441705\",\"GSM1441706\",\"GSM1441707\",\"GSM1441708\",\"GSM1441709\",\"GSM1441710\",\"GSM1441711\",\"GSM1441712\",\"GSM1441713\",\"GSM1441714\",\"GSM1441715\",\"GSM1441716\",\"GSM1441717\",\"GSM1441718\",\"GSM1441719\",\"GSM1441720\",\"GSM1441721\",\"GSM1441722\",\"GSM1441723\",\"GSM1441724\",\"GSM1441725\",\"GSM1441726\",\"GSM1441727\",\"GSM1441728\",\"GSM1441729\",\"GSM1441730\",\"GSM1441731\",\"GSM1441732\",\"GSM1441733\",\"GSM1441734\",\"GSM1441735\",\"GSM1441736\",\"GSM1441737\",\"GSM1441738\",\"GSM1441739\",\"GSM1441740\",\"GSM1441741\",\"GSM1441742\",\"GSM1441743\",\"GSM1441744\",\"GSM1441745\",\"GSM1441746\",\"GSM1441747\",\"GSM1441748\",\"GSM1441749\",\"GSM1441750\",\"GSM1441751\",\"GSM1441752\",\"GSM1441753\",\"GSM1441754\",\"GSM1441755\",\"GSM1441756\",\"GSM1441757\",\"GSM1441758\",\"GSM1441759\",\"GSM1441760\",\"GSM1441761\",\"GSM1441762\",\"GSM1441763\",\"GSM1441764\",\"GSM1441765\",\"GSM1441766\",\"GSM1441767\",\"GSM1441768\",\"GSM1441769\",\"GSM1441770\",\"GSM1441771\",\"GSM1441772\",\"GSM1441773\",\"GSM1441774\",\"GSM1441775\",\"GSM1441776\",\"GSM1441777\",\"GSM1441778\",\"GSM1441779\",\"GSM1441780\",\"GSM1441781\",\"GSM1441782\",\"GSM1441783\",\"GSM1441784\",\"GSM1441785\",\"GSM1441786\",\"GSM1441787\",\"GSM1441788\",\"GSM1441789\",\"GSM1441790\",\"GSM1441791\",\"GSM1441792\",\"GSM1441793\",\"GSM1441794\",\"GSM1441795\",\"GSM1441796\",\"GSM1441797\",\"GSM1441798\",\"GSM1441799\",\"GSM1441800\",\"GSM1441801\",\"GSM1441802\",\"GSM1441803\",\"GSM1441804\",\"GSM1441805\",\"GSM1441806\",\"GSM1441807\",\"GSM1441808\",\"GSM1441809\",\"GSM1441810\",\"GSM1441811\",\"GSM1441812\",\"GSM1441813\",\"GSM1441814\",\"GSM1441815\",\"GSM1441816\",\"GSM1441817\",\"GSM1441818\",\"GSM1441819\",\"GSM1441820\",\"GSM1441821\",\"GSM1441822\",\"GSM1441823\",\"GSM1441824\",\"GSM1441825\",\"GSM1441826\",\"GSM1441827\",\"GSM1441828\",\"GSM1441829\",\"GSM1445586\",\"GSM1445587\",\"GSM1445588\",\"GSM1445589\",\"GSM1445590\",\"GSM1445591\",\"GSM1445592\",\"GSM1445593\",\"GSM1445594\",\"GSM1445595\",\"GSM1445596\",\"GSM1445597\",\"GSM1445598\",\"GSM1445599\",\"GSM1445600\",\"GSM1445601\",\"GSM1445602\",\"GSM1445603\",\"GSM1445604\",\"GSM1445605\",\"GSM1445606\",\"GSM1445607\",\"GSM1445608\",\"GSM1445609\",\"GSM1445610\",\"GSM1445611\",\"GSM1445612\",\"GSM1445613\",\"GSM1445614\",\"GSM1445615\",\"GSM1445616\",\"GSM1445617\",\"GSM1445618\",\"GSM1445619\",\"GSM1445620\",\"GSM1445621\",\"GSM1445622\",\"GSM1445623\",\"GSM1445624\",\"GSM1445625\",\"GSM1445626\",\"GSM1445627\",\"GSM1445628\",\"GSM1445629\",\"GSM1445630\",\"GSM1445631\",\"GSM1445632\",\"GSM1445633\",\"GSM1445634\",\"GSM1445635\",\"GSM1445636\",\"GSM1445637\",\"GSM1445638\",\"GSM1445639\",\"GSM1445640\",\"GSM1445641\",\"GSM1445642\",\"GSM1445643\",\"GSM1445644\",\"GSM1445645\",\"GSM1445646\",\"GSM1445647\",\"GSM1445648\",\"GSM1445649\",\"GSM1445650\",\"GSM1445651\",\"GSM1445652\",\"GSM1445653\",\"GSM1445654\",\"GSM1445655\",\"GSM1445656\",\"GSM1445657\",\"GSM1445658\",\"GSM1445659\",\"GSM1445660\",\"GSM1445661\",\"GSM1445662\",\"GSM1445663\",\"GSM1445664\",\"GSM1445665\",\"GSM1445666\",\"GSM1445667\",\"GSM1445668\",\"GSM1445669\",\"GSM1445670\",\"GSM1445671\",\"GSM1445672\",\"GSM1445673\",\"GSM1445674\",\"GSM1445675\",\"GSM1445676\",\"GSM1445677\",\"GSM1445678\",\"GSM1445679\",\"GSM1445680\",\"GSM1445681\",\"GSM1445682\",\"GSM1445683\",\"GSM1445684\",\"GSM1445685\",\"GSM1445686\",\"GSM1445687\",\"GSM1445688\",\"GSM1445689\",\"GSM1445690\",\"GSM1445691\",\"GSM1445692\",\"GSM1445693\",\"GSM1445694\",\"GSM1445695\",\"GSM1445696\",\"GSM1445697\",\"GSM1445698\",\"GSM1445699\",\"GSM1445700\",\"GSM1445701\",\"GSM1445702\",\"GSM1445703\",\"GSM1445704\",\"GSM1445705\",\"GSM1445706\",\"GSM1445707\",\"GSM1445708\",\"GSM1445709\",\"GSM1445710\",\"GSM1445711\",\"GSM1445712\",\"GSM1445713\",\"GSM1445714\",\"GSM1445715\",\"GSM1445716\",\"GSM1445717\",\"GSM1445718\",\"GSM1445719\",\"GSM1445720\",\"GSM1445721\",\"GSM1445722\",\"GSM1445723\",\"GSM1445724\",\"GSM1445725\",\"GSM1445726\",\"GSM1445727\",\"GSM1445728\",\"GSM1445729\",\"GSM1445730\",\"GSM1445731\",\"GSM1445732\",\"GSM1445733\",\"GSM1445734\",\"GSM1445735\",\"GSM1445736\",\"GSM1445737\",\"GSM1445738\",\"GSM1445739\",\"GSM1445740\",\"GSM1445741\",\"GSM1445742\",\"GSM1445743\",\"GSM1445744\",\"GSM1445745\",\"GSM1445746\",\"GSM1445747\",\"GSM1445748\",\"GSM1445749\",\"GSM1445750\",\"GSM1445751\",\"GSM1445752\",\"GSM1445753\",\"GSM1445754\",\"GSM1445755\",\"GSM1445756\",\"GSM1445757\",\"GSM1445758\",\"GSM1445759\",\"GSM1445760\",\"GSM1445761\",\"GSM1445762\",\"GSM1445763\",\"GSM1445764\",\"GSM1445765\",\"GSM1445766\",\"GSM1445767\",\"GSM1445768\",\"GSM1445769\",\"GSM1445770\",\"GSM1445771\",\"GSM1445772\",\"GSM1445773\",\"GSM1445774\",\"GSM1445775\",\"GSM1445776\",\"GSM1445777\",\"GSM1445778\",\"GSM1445779\",\"GSM1445780\",\"GSM1445781\",\"GSM1445782\",\"GSM1445783\",\"GSM1445784\",\"GSM1445785\",\"GSM1445786\",\"GSM1445787\",\"GSM1445788\",\"GSM1445789\",\"GSM1445790\",\"GSM1445791\",\"GSM1445792\",\"GSM1445793\",\"GSM1445794\",\"GSM1445795\",\"GSM1445796\",\"GSM1445797\",\"GSM1445798\",\"GSM1445799\",\"GSM1445800\",\"GSM1445801\",\"GSM1445802\",\"GSM1445803\",\"GSM1445804\",\"GSM1445805\",\"GSM1445806\",\"GSM1445807\",\"GSM1445808\",\"GSM1445809\",\"GSM1445810\",\"GSM1445811\",\"GSM1445812\",\"GSM1445813\",\"GSM1445814\",\"GSM1445815\",\"GSM1445816\",\"GSM1445817\",\"GSM1445818\",\"GSM1445819\",\"GSM1445820\",\"GSM1445821\"],\n"
    			+ "  \"studyCohort\": [\"age_group\", \"gender\"],\n"
    			+ "  \"studyVariables\": {\n"
    			+ "    \"vaccine\": [\"Fluvirin\"],\n"
    			+ "    \"study\": [\"SDY400\", \"SDY404\", \"SDY520\"],\n"
    			+ "    \"age_group\": [\"old\", \"young\"]\n"
    			+ "  },\n"
    			+ "  \"modelTime\": true,\n"
    			+ "  \"analysisGroups\": {\n"
    			+ "    \"group1\": [0],\n"
    			+ "    \"group2\": [2,  3,  4,  5,  7,  8,  9, 24, 28, 32, 35, 36, 43]\n"
    			+ "  },\n"
    			+ "  \"platformCorrection\": true, \n"
    			+ "  \"variableGenes\": true\n"
    			+ "}");
    	outputJSON(rtn);
    }
    
    @Test
    public void testPathwayAnalysisForGenes() throws Exception{
    	String url = HOST_URL + "/analysis/pathways";
    	List<String> genes = new ArrayList<>(Arrays.asList("LRRN3", "SERPINE2", "ACTA2", "CCR7", "SYT11", "ADRB2"));
    	ObjectMapper mapper = new ObjectMapper();
    	String json = mapper.writeValueAsString(genes);
    	String rtn = callHttp(url, HTTP_POST, json);
    	outputJSON(rtn);
    }
    
    @Test
    public void testGetBiosampleMetadata() throws Exception {
    	String url = HOST_URL + "/metadata/biosamples";
    	String rtn = callHttp(url, HTTP_GET, null);
    	System.out.println(rtn);
    }
    
    private String outputJSON(String json) throws JsonProcessingException, IOException {
        ObjectMapper mapper = new ObjectMapper();
        Object obj = mapper.readValue(json, Object.class);
        String rtn = mapper.writerWithDefaultPrettyPrinter().writeValueAsString(obj);
        System.out.println(rtn);
        return rtn;
    }

    protected String callHttp(String url,
                              String type,
                              String query) throws IOException {
        HttpMethod method = null;
        HttpClient client = null;
        if (type.equals(HTTP_POST)) {
            method = new PostMethod(url);
            client = initializeHTTPClient((PostMethod) method, query);
        } else {
            method = new GetMethod(url); // Default
            client = new HttpClient();
        }
        method.setRequestHeader("Accept", "application/json");
        method.setRequestHeader("content-type", "application/json");
        int responseCode = client.executeMethod(method);
        if (responseCode == HttpStatus.SC_OK) {
            InputStream is = method.getResponseBodyAsStream();
            return readMethodReturn(is);
        } else {
            System.err.println("Error from server: " + method.getResponseBodyAsString());
            System.out.println("Response code: " + responseCode);
            throw new IllegalStateException(method.getResponseBodyAsString());
        }
    }

    protected String readMethodReturn(InputStream is) throws IOException {
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader reader = new BufferedReader(isr);
        StringBuilder builder = new StringBuilder();
        String line = null;
        while ((line = reader.readLine()) != null)
            builder.append(line).append("\n");
        reader.close();
        isr.close();
        is.close();
        // Remove the last new line
        String rtn = builder.toString();
        // Just in case an empty string is returned
        if (rtn.length() == 0)
            return rtn;
        return rtn.substring(0, rtn.length() - 1);
    }

    private HttpClient initializeHTTPClient(PostMethod post, String query) throws UnsupportedEncodingException {
        RequestEntity entity = new StringRequestEntity(query, "text/plain", "UTF-8");
        //        RequestEntity entity = new StringRequestEntity(query, "application/XML", "UTF-8");
        post.setRequestEntity(entity);
        //        post.setRequestHeader("Accept", "application/JSON, application/XML, text/plain");
        post.setRequestHeader("Accept", "application/json");
        HttpClient client = new HttpClient();
        return client;
    }

}
