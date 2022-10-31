import os, sys
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import tensorflow as tf
tf.get_logger().setLevel('ERROR')
import numpy as np
from numpy import *
##############################################################################

##############################################################################
## System Paths ##
#path                 = './other_soft/DeepPE_modified/DeepPE/DeepPE_test_input.txt'

## Run Parameters ##
TEST_NUM_SET         = [0] #to be ignored - (leaving the original code structure as it is)
PACKAGE_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
MODEL_DIR = PACKAGE_DIRECTORY + '/DeepPE_weights/'
best_model_path_list = [MODEL_DIR]
#

#output path
#testbook = xlsxwriter.Workbook('./other_soft/DeepPE_modified/DeepPE/seri7_DeepPE_example_output.xlsx')

length = 47

# Model
class Model(object):
    def __init__(self, filter_size, filter_num, length, node_1 = 80, node_2 = 60, l_rate = 0.005, bio_num = 20):
        L_filter_num = 4
        bio_num = 20
        self.inputs         = tf.compat.v1.placeholder(tf.float32, [None, 1, length, L_filter_num])
        self.mod_inputs     = tf.compat.v1.placeholder(tf.float32, [None, 1, length, L_filter_num])
        self.bios           = tf.compat.v1.placeholder(tf.float32, [None, bio_num])
        self.targets        = tf.compat.v1.placeholder(tf.float32, [None, 1])
        self.is_training    = tf.compat.v1.placeholder(tf.bool)
        
        def create_new_conv_layer(input_data, num_input_channels, num_filters, filter_shape, stride_shape, name):
            conv_filt_shape = [filter_shape[0], filter_shape[1], num_input_channels, num_filters]
            weights   = tf.Variable(tf.compat.v1.truncated_normal(conv_filt_shape, stddev=0.03),name=name+'_W')
            bias      = tf.Variable(tf.compat.v1.truncated_normal([num_filters]), name=name+'_b')
            out_layer = tf.nn.conv2d(input_data, weights, [1, stride_shape[0], stride_shape[1], 1], padding='VALID')
            out_layer += bias
            out_layer = tf.compat.v1.layers.dropout(tf.nn.relu(out_layer), 0.3, self.is_training)
            return out_layer
            
        stride = 1
        L_pe = create_new_conv_layer(self.inputs, L_filter_num, filter_num[0], [1, filter_size[0]], [1, stride], name='pe')
        L_mod_pe = create_new_conv_layer(self.mod_inputs, L_filter_num, filter_num[0], [1, filter_size[0]], [1, stride], name='pe_mod')
        
        layer_node = int((length-filter_size[0])/stride)+1
        node_num_0 = layer_node*filter_num[0]
        
        L_flatten_pe  = tf.reshape(L_pe, [-1, node_num_0])
        L_flatten_mod_pe  = tf.reshape(L_mod_pe, [-1, node_num_0])
        L_flatten = tf.concat([L_flatten_pe, L_flatten_mod_pe, self.bios], 1, name='concat')
        
        with tf.compat.v1.variable_scope('Fully_Connected_Layer1'):
            W_fcl1       = tf.compat.v1.get_variable("W_fcl1", shape=[node_num_0*2+bio_num, node_1])
            B_fcl1       = tf.compat.v1.get_variable("B_fcl1", shape=[node_1])
            L_fcl1_pre   = tf.nn.bias_add(tf.matmul(L_flatten, W_fcl1), B_fcl1)
            L_fcl1       = tf.nn.relu(L_fcl1_pre)
            L_fcl1_drop  = tf.compat.v1.layers.dropout(L_fcl1, 0.3, self.is_training)
            
        if node_2 != 0:
            with tf.compat.v1.variable_scope('Fully_Connected_Layer2'):
                W_fcl2       = tf.compat.v1.get_variable("W_fcl2", shape=[node_1, node_2])
                B_fcl2       = tf.compat.v1.get_variable("B_fcl2", shape=[node_2])
                L_fcl2_pre   = tf.nn.bias_add(tf.matmul(L_fcl1_drop, W_fcl2), B_fcl2)
                L_fcl2       = tf.nn.relu(L_fcl2_pre)
                L_fcl2_drop  = tf.layers.dropout(L_fcl2, 0.3, self.is_training)
            with tf.variable_scope('Output_Layer'):
                W_out        = tf.get_variable("W_out", shape=[node_2, 1])
                B_out        = tf.get_variable("B_out", shape=[1])
                self.outputs = tf.nn.bias_add(tf.matmul(L_fcl2_drop, W_out), B_out)        
        else:
            with tf.compat.v1.variable_scope('Output_Layer'):
                W_out        = tf.compat.v1.get_variable("W_out", shape=[node_1, 1])
                B_out        = tf.compat.v1.get_variable("B_out", shape=[1])
                self.outputs = tf.nn.bias_add(tf.matmul(L_fcl1_drop, W_out), B_out)       
        self.obj_loss    = tf.reduce_mean(tf.square(self.targets - self.outputs))
        self.optimizer   = tf.compat.v1.train.AdamOptimizer(l_rate).minimize(self.obj_loss)

def preprocess_seq(data):
    #print("Start preprocessing the sequence done 2d")
    global length
    DATA_X = np.zeros((len(data),1,length,4), dtype=int)
    #print(np.shape(data), len(data), length)
    for l in range(len(data)):
        for i in range(length):
            try: data[l][i]
            except: print(data[l], i, length, len(data))
            if data[l][i]in "Aa":    DATA_X[l, 0, i, 0] = 1
            elif data[l][i] in "Cc": DATA_X[l, 0, i, 1] = 1
            elif data[l][i] in "Gg": DATA_X[l, 0, i, 2] = 1
            elif data[l][i] in "Tt": DATA_X[l, 0, i, 3] = 1
            elif data[l][i] in "Xx": pass
            else:
                print("Non-ATGC character " + data[l])
                print(i)
                print(data[l][i])
                sys.exit()
        #loop end: i
    #loop end: l
    #print("Preprocessing the sequence done")
    return DATA_X
#def end: preprocess_seq


def getfile_inference(wide_arg, edit_arg, pbs_len_arg, rtt_len_arg, pbs_rtt_len_arg, Tm1_arg, Tm2_arg, Tm3_arg,Tm4_arg,deltaTm_arg,GCcount1_arg,GCcount2_arg,GCcount3_arg,GCcontent1_arg, GCcontent2_arg,GCcontent3_arg,mfe_1_arg, mfe_2_arg,mfe_3_arg, mfe_4_arg,mfe_5_arg, deep_sp_cas9_arg):
    seq     = [wide_arg]
    mod_seq     = [edit_arg]
    bio = []
    tmp_bio = [pbs_len_arg, rtt_len_arg, pbs_rtt_len_arg, Tm1_arg, Tm2_arg, Tm3_arg,Tm4_arg,deltaTm_arg,GCcount1_arg,GCcount2_arg,GCcount3_arg,GCcontent1_arg, GCcontent2_arg,GCcontent3_arg,mfe_1_arg, mfe_2_arg,mfe_3_arg, mfe_4_arg,mfe_5_arg, deep_sp_cas9_arg]
    bio.append(tmp_bio)
    processed_seq = preprocess_seq(seq)
    processed_mod_seq = preprocess_seq(mod_seq)
    return processed_seq, seq, processed_mod_seq, mod_seq, bio
#def end: getseq

# Test Model
def Model_Inference(sess, TEST_X, model, TEST_mod_X=None, TEST_bio=None):
    test_batch = 500
    TEST_Z = zeros((TEST_X.shape[0], 1), dtype=float)
    for i in range(int(ceil(float(TEST_X.shape[0])/float(test_batch)))):
        Dict = {model.inputs: TEST_X[i*test_batch:(i+1)*test_batch],
                    model.mod_inputs: TEST_mod_X[i*test_batch:(i+1)*test_batch],
                    model.bios: TEST_bio[i*test_batch:(i+1)*test_batch],
                    model.is_training: False}
        TEST_Z[i*test_batch:(i+1)*test_batch] = sess.run([model.outputs], feed_dict=Dict)[0]
    #testval_row = 0
    #testval_col = 3
    #sheet_index = 0
    #for test_value in (TEST_Z):
    #    testvalsheet[sheet_index].write(testval_row, testval_col, test_value[0])
    #    testval_row += 1
    return TEST_Z[0]

def main(wide_arg, edit_arg, pbs_len_arg, rtt_len_arg, pbs_rtt_len_arg, Tm1_arg, Tm2_arg, Tm3_arg,Tm4_arg,deltaTm_arg,GCcount1_arg,GCcount2_arg,GCcount3_arg,GCcontent1_arg, GCcontent2_arg,GCcontent3_arg,mfe_1_arg, mfe_2_arg,mfe_3_arg, mfe_4_arg,mfe_5_arg, deep_sp_cas9_arg):
    #TensorFlow config
    conf                                = tf.compat.v1.ConfigProto()
    conf.gpu_options.allow_growth       = True
    os.environ['CUDA_VISIBLE_DEVICES']  = '0'
    best_model_cv                       = 0.0
    best_model_list                     = []

    TEST_X = []
    TEST_mod_X = []
    TEST_bio = []
    testsheet = []
    for TEST_NUM_index in range(len(TEST_NUM_SET)):
        TEST_NUM = TEST_NUM_SET[TEST_NUM_index]
        #testsheet.append([testbook.add_worksheet('{}'.format(TEST_NUM))])
        tmp_X, pre_X, tmp_mod_X, pre_mod_X, bio = getfile_inference(wide_arg, edit_arg, pbs_len_arg, rtt_len_arg, pbs_rtt_len_arg, Tm1_arg, Tm2_arg, Tm3_arg,Tm4_arg,deltaTm_arg,GCcount1_arg,GCcount2_arg,GCcount3_arg,GCcontent1_arg, GCcontent2_arg,GCcontent3_arg,mfe_1_arg, mfe_2_arg,mfe_3_arg, mfe_4_arg,mfe_5_arg, deep_sp_cas9_arg)
        TEST_X.append(tmp_X)
        TEST_mod_X.append(tmp_mod_X)
        TEST_bio.append(bio)
        #test_row = 0
        #for index_X in range(np.shape(pre_X)[0]):
        #    testsheet[-1][-1].write(test_row, 0, pre_X[index_X])
        #    testsheet[-1][-1].write(test_row, 1, pre_mod_X[index_X])
        #    testsheet[-1][-1].write(test_row, 2, repr(bio[index_X]))
        #    test_row += 1

    for best_model_path in best_model_path_list:
        for modelname in os.listdir(best_model_path):
            if "meta" in modelname:
                best_model_list.append(modelname[:-5])

    for index in range(len(best_model_list)):
        best_model_path = best_model_path_list[index]
        best_model      = best_model_list[index]
        valuelist       = best_model.split('-')
        fulllist        = []
        for value in valuelist:
            try:
                value=int(value)
            except:
                try:    value=float(value)
                except: pass
            fulllist.append(value)
        #print(fulllist[2:])
        filter_size_1, filter_size_2, filter_size_3, filter_num_1, filter_num_2, filter_num_3, l_rate, load_episode, node_1 = fulllist[2:]
        filter_size = [filter_size_1, filter_size_2, filter_size_3]
        filter_num  = [filter_num_1, filter_num_2, filter_num_3]
        args = [filter_size, filter_num, l_rate, 0, None, node_1]
        # Loading the model with the best validation score and test
        tf.compat.v1.reset_default_graph()
        with tf.compat.v1.Session(config=conf) as sess:
            sess.run(tf.compat.v1.global_variables_initializer())
            model = Model(filter_size, filter_num, length, node_1, 0, args[2])
            saver = tf.compat.v1.train.Saver()
            saver.restore(sess, best_model_path+"/PreTrain-Final-{}-{}-{}-{}-{}-{}-{}-{}-{}".format(args[0][0], args[0][1], args[0][2], args[1][0], args[1][1], args[1][2], args[2], load_episode, args[5]))
            result = Model_Inference(sess, TEST_X[0], model, TEST_mod_X=TEST_mod_X[0], TEST_bio=TEST_bio[0])[0]
            #testbook.close()
            return result


if __name__ == '__main__':
    main()